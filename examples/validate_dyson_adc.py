"""Molecular smoke/comparison runs, not an external ADC(3) benchmark."""
import json
from pathlib import Path
import pyscf
from pyscf import gto, scf
from nondiagonal_ept import EPT

cases = [
    ('HF', 'F 0 0 0; H 0 0 .9168', 'cc-pvdz', 1, [4, 2], [5]),
    ('HF', 'F 0 0 0; H 0 0 .9168', 'cc-pvtz', 1, [4, 2], [5]),
    ('H2O', 'O 0 0 0; H 0 -.7586 .5043; H 0 .7586 .5043', 'cc-pvdz', 1, [4, 3], [5]),
    ('N2', 'N 0 0 0; N 0 0 1.1136', 'cc-pvdz', 2, [6, 4, 3], [7]),
]
report = {'pyscf_version': pyscf.__version__, 'cases': []}
for name, atom, basis, frozen, ip, ea in cases:
    mf = scf.RHF(gto.M(atom=atom, unit='Angstrom', basis=basis, verbose=0)).run(conv_tol=1e-12)
    # Dyson Hamiltonians are the same in both sectors: reuse each for both.
    for method in ('ADC(2)', '3+', 'ADC(3)'):
        calc = EPT(mf, method, frozen=frozen)
        data = dict(molecule=name, atom=atom, unit='Angstrom', basis=basis,
                    frozen=frozen, method=method, hf_energy_hartree=mf.e_tot, poles=[])
        for sector, targets in [('ip', ip), ('ea', ea)]:
            calc.hamiltonian.sector = sector
            for pole in calc.kernel(targets):
                data['poles'].append(dict(sector=sector, target=pole.target,
                    omega_hartree=pole.energy, binding_energy_ev=pole.binding_energy_ev,
                    PS=pole.strength, residual=pole.residual))
        if hasattr(calc.hamiltonian,'static_diagnostics'):
            data['static_diagnostics'] = {k:v for k,v in calc.hamiltonian.static_diagnostics.items()
                                         if k not in ('dynamic_density','correlation_density')}
        report['cases'].append(data)
        print(name, basis, method, 'complete', flush=True)
Path('docs/dyson_adc_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
