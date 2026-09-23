"""Conventional Dyson ADC(2), ADC(3)-DEM and strict third-order 3+."""
from pyscf import gto, scf
from nondiagonal_ept import run_methods

mol = gto.M(
    atom="""
    F   0.0000   0.0000   0.0000
    H   0.0000   0.0000   0.9168
    """,
    unit='Angstrom', basis='cc-pvdz', verbose=0,
)
mf = scf.RHF(mol).run(conv_tol=1e-12)
for sector, targets in [('ip', [4, 2]), ('ea', [5])]:
    results = run_methods(mf, ['ADC(2)', 'ADC(3)', '3+'], frozen=1,
                          sector=sector, targets=targets, tol=1e-9)
    for name, poles in results.items():
        label = 'ADC(2) / ND2' if name == 'ND2' else name
        for p in poles:
            print(f'{label:14s} {sector.upper()} MO {p.target}: '
                  f'{p.binding_energy_ev:12.6f} eV, PS = {p.strength:.6f}, '
                  f'residual = {p.residual:.2e}')
