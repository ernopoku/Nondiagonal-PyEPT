"""Measure complete EPT setup and pole solves, excluding the shared RHF.

Run this script with both old and updated installations for a timing comparison:
  python examples/benchmark_non_dyson.py --basis cc-pvtz --output timing.json
Thread settings are applied before importing numerical libraries.
"""
import argparse
import json
import os
import platform
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--basis', default='cc-pvtz')
    parser.add_argument('--molecule', choices=['HF', 'H2O', 'N2'], default='HF')
    parser.add_argument('--sector', choices=['ip', 'ea'], default='ip')
    parser.add_argument('--method', choices=['NRL3', 'nD-NRL3'], default='nD-NRL3')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--memory-mb', type=float, default=16000)
    parser.add_argument('--verbose', type=int, default=0)
    parser.add_argument('--output', type=Path, default=Path('non_dyson_timing.json'))
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'VECLIB_MAXIMUM_THREADS'):
        os.environ[key] = str(args.threads)
    import numpy as np
    import scipy
    import pyscf
    from pyscf import gto, scf, lib
    import nondiagonal_ept
    from nondiagonal_ept import EPT
    lib.num_threads(args.threads)
    cases = {
        'HF': ('F 0 0 0; H 0 0 .9168', 1, [4, 2], [5]),
        'H2O': ('O 0 0 0; H 0 -.7586 .5043; H 0 .7586 .5043', 1, [4, 3], [5]),
        'N2': ('N 0 0 0; N 0 0 1.1136', 2, [6, 4, 3], [7]),
    }
    atom, frozen, ip, ea = cases[args.molecule]
    mf = scf.RHF(gto.M(atom=atom, basis=args.basis, unit='Angstrom',
                      verbose=args.verbose, max_memory=args.memory_mb))
    mf.conv_tol = 1e-12
    mf.kernel()
    if not mf.converged:
        raise RuntimeError('RHF failed to converge')
    wall, cpu = time.perf_counter(), time.process_time()
    calculation = EPT(mf, args.method, frozen=frozen, sector=args.sector,
                      max_memory_mb=args.memory_mb)
    setup_wall, setup_cpu = time.perf_counter()-wall, time.process_time()-cpu
    wall, cpu = time.perf_counter(), time.process_time()
    poles = calculation.kernel(ip if args.sector == 'ip' else ea, tol=1e-9)
    solve_wall, solve_cpu = time.perf_counter()-wall, time.process_time()-cpu
    h = calculation.hamiltonian
    result = dict(
        molecule=args.molecule, basis=args.basis, atom=atom, frozen=frozen,
        sector=args.sector, method=calculation.method, threads=args.threads,
        hf_energy=float(mf.e_tot), dimension=h.shape[0],
        setup_seconds=setup_wall, solve_seconds=solve_wall,
        total_seconds=setup_wall+solve_wall, setup_cpu_seconds=setup_cpu,
        solve_cpu_seconds=solve_cpu,
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         scipy=scipy.__version__, pyscf=pyscf.__version__,
                         platform=platform.platform(), processor=platform.processor(),
                         source=str(Path(nondiagonal_ept.__file__).resolve())),
        static_diagnostics=getattr(h, 'static_diagnostics', {}),
        poles=[dict(target=p.target, energy=p.energy,
                    binding_energy_ev=p.binding_energy_ev, PS=p.strength,
                    residual=p.residual) for p in poles],
    )
    if hasattr(h, 'static_opposite'):
        # Small matrix included for full-correction comparison, not just roots.
        result['static_opposite'] = h.static_opposite.tolist()
        result['static_residuals'] = h.static_solve_residuals.tolist()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(f'{args.molecule}/{args.basis} {args.method} {args.sector}: '
          f'setup {setup_wall:.3f} s, poles {solve_wall:.3f} s; {args.output}')


if __name__ == '__main__':
    main()
