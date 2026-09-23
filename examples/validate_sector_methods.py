"""Reproducible small-basis comparisons; not an external NRL3-ISR benchmark."""
import json
from pathlib import Path
import numpy as np
import pyscf
from pyscf import gto, scf
from nondiagonal_ept import EPT, SectorEPT

cases = [
    ('HF', 'F 0 0 0; H 0 0 .9168', '6-31g', 1),
    ('H2O', 'O 0 0 0; H 0 -.7586 .5043; H 0 .7586 .5043', '6-31g', 1),
    ('N2', 'N 0 0 0; N 0 0 1.1136', 'sto-3g', 2),
]
report = {'pyscf_version': pyscf.__version__, 'cases': []}
for name, atom, basis, frozen in cases:
    mf = scf.RHF(gto.M(atom=atom, unit='Angstrom', basis=basis, verbose=0)).run(conv_tol=1e-12)
    dyson = EPT(mf, 'NRL3', frozen=frozen).dense_spectrum(max_dimension=1200)
    # Separate parent roots by the gap, only for these stable closed-shell cases.
    gap_center = .5*(mf.mo_energy[mf.mo_occ == 2][-1]+mf.mo_energy[mf.mo_occ == 0][0])
    for sector in ('ip', 'ea'):
        parent = sorted([p for p in dyson if (p.energy < gap_center) == (sector == 'ip')],
                        key=lambda p: p.energy, reverse=sector == 'ip')[:3]
        row = dict(molecule=name, atom=atom, basis=basis, unit='Angstrom', frozen=frozen, sector=sector)
        row['NRL3'] = [{'binding_energy_ev': p.binding_energy_ev, 'PS': p.strength} for p in parent]
        for method in ('nD-ADC(3)', 'NRL3-ISR(3)'):
            calc = SectorEPT(mf, method, frozen=frozen, sector=sector)
            poles = calc.kernel(nroots=3)
            row[method] = [{'binding_energy_ev': p.binding_energy_ev, 'PS': p.strength,
                            'residual': p.residual} for p in poles]
            if hasattr(calc, 'diagnostics'):
                row['diagnostics'] = calc.diagnostics
        report['cases'].append(row)
path = Path('docs/sector_comparison.json')
path.write_text(json.dumps(report, indent=2, default=lambda v: v.item() if isinstance(v, np.generic) else v)+'\n')
print(f'Wrote {path}')
