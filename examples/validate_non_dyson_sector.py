"""Validate sector-projected nD-NRL3 vs NRL3 on supplied neutral geometries.

Each molecule/basis runs in a subprocess with a timeout and checkpointed output.
Results are numerical comparisons to NRL3, not experimental EA/IP benchmarks.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

# User-supplied coordinates, in Angstrom. Neutral singlet RHF; frozen 1s cores.
CASES = {
    'LiH': ('Li 0 0 .401144; H 0 0 -1.203431', 1),
    'Li2': ('Li 0 0 1.374955; Li 0 0 -1.374955', 2),
    'LiF': ('F 0 0 .398667; Li 0 0 -1.196001', 2),
    'BeO': ('Be 0 0 -.904151; O 0 0 .452076', 2),
    'C3': ('C 0 0 0; C 0 0 1.299192; C 0 0 -1.299192', 3),
}
BASES = ['cc-pvdz', 'aug-cc-pvdz', 'cc-pvtz', 'aug-cc-pvtz']


def worker(args):
    # Set thread controls before any scientific imports.
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'VECLIB_MAXIMUM_THREADS'):
        os.environ[key] = str(args.threads)
    import platform
    import numpy as np
    import scipy
    import pyscf
    from pyscf import gto, scf, lib
    from nondiagonal_ept.integrals import from_pyscf
    from nondiagonal_ept.blocks import Hamiltonian
    from nondiagonal_ept.non_dyson import StaticNRL3
    from nondiagonal_ept.solver import davidson, make_pole
    lib.num_threads(args.threads)
    molecule, basis = args.worker
    atom, frozen = CASES[molecule]
    record = dict(molecule=molecule, basis=basis, atom=atom, unit='Angstrom',
                  charge=0, spin=0, frozen=frozen, threads_requested=args.threads,
                  memory_mb=args.memory_mb, scf_tol=1e-12, pole_tol=1e-9,
                  static_tol=1e-10, formulation='sector-projected-v2', sectors=[],
                  environment=dict(python=platform.python_version(), numpy=np.__version__,
                                   scipy=scipy.__version__, pyscf=pyscf.__version__))
    import inspect
    record['non_dyson_source_sha256'] = hashlib.sha256(
        Path(inspect.getfile(StaticNRL3)).read_bytes()).hexdigest()
    path = args.output / f'{molecule}-{basis}.json'
    started = time.perf_counter()
    def save(stage):
        record['stage'] = stage
        record['elapsed_seconds'] = time.perf_counter() - started
        path.write_text(json.dumps(record, indent=2) + '\n')
        print(molecule, basis, stage, flush=True)
    try:
        save('SCF')
        mol = gto.M(atom=atom, basis=basis, charge=0, spin=0, unit='Angstrom',
                    verbose=0, max_memory=args.memory_mb)
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.max_cycle = 200
        mf.kernel()
        if not mf.converged:
            mf = mf.newton().run(conv_tol=1e-12, max_cycle=100)
            # Return ordinary RHF with the converged reference orbitals.
            mf = mf.undo_soscf()
        if not mf.converged:
            raise RuntimeError('RHF did not converge')
        save('RHF stability')
        repairs = 0
        while True:
            mo, _, stable_internal, stable_external = mf.stability(
                internal=True, external=True, return_status=True)
            if stable_internal or repairs == 3:
                break
            mf.kernel(dm0=mf.make_rdm1(mo, mf.mo_occ))
            if not mf.converged:
                raise RuntimeError('RHF stability restart did not converge')
            repairs += 1
        record['reference'] = dict(energy_hartree=float(mf.e_tot),
                                  stable_internal=bool(stable_internal),
                                  stable_external=bool(stable_external),
                                  internal_restarts=repairs)
        if not stable_internal:
            raise RuntimeError('Internally unstable RHF after three restarts')
        # External instability is reported, not silently replaced with UHF.
        save('integrals')
        ints = from_pyscf(mf, frozen=frozen, max_memory_mb=args.memory_mb)
        nocc = mol.nelectron // 2
        targets = {'ip': list(range(nocc-1, max(frozen-1,nocc-3), -1)),
                   'ea': [nocc, nocc+1]}
        for sector in ('ip', 'ea'):
            item = dict(sector=sector, targets=targets[sector], poles=[])
            record['sectors'].append(item)
            try:
                save(f'{sector}: Hamiltonians')
                parent = Hamiltonian(ints, 'NRL3', sector)
                nd = StaticNRL3(ints, sector, max_memory_mb=args.memory_mb,
                               static_space='sector', static_tol=1e-10)
                item['static_diagnostics'] = nd.static_diagnostics
                item['dimensions'] = {'NRL3': parent.shape[0], 'nD-NRL3': nd.shape[0]}
                seen = {'NRL3': [], 'nD-NRL3': []}
                for target in targets[sector]:
                    entry = dict(target=target)
                    item['poles'].append(entry)
                    try:
                        found = []
                        for label, ham in [('NRL3', parent), ('nD-NRL3', nd)]:
                            save(f'{sector}: MO {target} {label}')
                            original = ints.original_mos[ints.spatial[ham.simple]]
                            index = int(np.flatnonzero(original == target)[0])
                            e, x, residual, iterations = davidson(
                                ham, index, tol=1e-9, max_cycle=300, max_space=60)
                            for previous_target, previous_vector in seen[label]:
                                if abs(previous_vector @ x) > 1-1e-7:
                                    raise RuntimeError(f'{label} target {target} duplicates root for target {previous_target}')
                            seen[label].append((target, x.copy()))
                            pole = make_pole(ham, e, x, residual, iterations, target)
                            entry[label] = dict(binding_energy_eV=pole.binding_energy_ev,
                                                PS=pole.strength, residual=residual,
                                                iterations=iterations)
                            found.append(pole)
                        a, b = found
                        entry['difference_eV'] = b.binding_energy_ev-a.binding_energy_ev
                        entry['dyson_overlap_squared'] = float((a.normalized_dyson_mo@b.normalized_dyson_mo)**2)
                    except Exception as exc:
                        entry['error'] = repr(exc)
                del parent, nd
            except Exception as exc:
                item['error'] = repr(exc)
            save(f'{sector}: complete')
        save('complete')
    except Exception as exc:
        record['error'] = repr(exc)
        save('failed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--molecules', nargs='+', choices=list(CASES), default=list(CASES))
    parser.add_argument('--bases', nargs='+', default=BASES)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--memory-mb', type=float, default=16000)
    parser.add_argument('--timeout', type=float, default=900, help='Seconds per molecule/basis')
    parser.add_argument('--output', type=Path, default=Path('nd-sector-validation'))
    parser.add_argument('--worker', nargs=2, help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.worker:
        worker(args)
        return
    summary = []
    for molecule in args.molecules:
        for basis in args.bases:
            name = f'{molecule}-{basis}'
            print('START', name, flush=True)
            command = [sys.executable, str(Path(__file__).resolve()), '--worker', molecule, basis,
                       '--threads', str(args.threads), '--memory-mb', str(args.memory_mb),
                       '--output', str(args.output.resolve())]
            status = dict(molecule=molecule, basis=basis)
            with (args.output/f'{name}.log').open('w') as log:
                try:
                    process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                             timeout=args.timeout, check=False)
                    status['exit_code'] = process.returncode
                except subprocess.TimeoutExpired:
                    status['timeout_seconds'] = args.timeout
            result_path = args.output/f'{name}.json'
            if result_path.exists():
                status['result'] = json.loads(result_path.read_text())
            result = status.get('result', {})
            status['validation_ok'] = (status.get('exit_code') == 0
                and result.get('stage') == 'complete' and 'error' not in result
                and all('error' not in sector and all('error' not in pole for pole in sector['poles'])
                        for sector in result.get('sectors', [])))
            summary.append(status)
            (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
            print('DONE' if status['validation_ok'] else 'FAILED', name,
                  {k:v for k,v in status.items() if k != 'result'}, flush=True)
    if any(not case['validation_ok'] for case in summary):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
