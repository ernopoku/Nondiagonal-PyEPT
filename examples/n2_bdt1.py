"""N2/cc-pVDZ BD-T1 with two fixed RHF core orbitals."""
import json
import numpy as np
from pyscf import gto, scf
from nondiagonal_ept import EPT

mol = gto.M(atom='N 0 0 0; N 0 0 1.1136', unit='Angstrom',
            basis='cc-pvdz', verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)
calculation = EPT(mf, 'BD-T1', frozen=2, brueckner_tol=1e-8)
reference = [-0.6139503411, -0.5491151633, -0.6742780381]
strengths = [0.9292954916, 0.9059496833, 0.8545204419]
report = {'molecule': 'N2', 'bond_angstrom': 1.1136, 'basis': 'cc-pvdz',
          'frozen': 2, 'core_convention': 'fixed input RHF orbitals',
          'bd_energy_hartree': calculation.brueckner.e_tot,
          'singles_norm': float(np.linalg.norm(calculation.brueckner.t1)),
          'core_coefficient_max_change': float(np.max(np.abs(
              calculation.brueckner.mo_coeff[:,:2] - mf.mo_coeff[:,:2]))),
          'poles': []}
for p, ref, z in zip(calculation.kernel([6,4,3]), reference, strengths):
    report['poles'].append({'original_mo': p.target, 'energy': p.energy,
        'reference_energy': ref, 'difference_hartree': p.energy-ref,
        'strength': p.strength, 'reference_strength': z, 'residual': p.residual})
print(json.dumps(report, indent=2))
