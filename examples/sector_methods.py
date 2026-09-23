"""Separate IP/EA methods on one RHF reference; small basis for the dense prototype.

Run: python examples/sector_methods.py
Use nD-ADC(3) alone for larger bases; NRL3-ISR(3) is experimental and dense.
"""
from pyscf import gto, scf
from nondiagonal_ept import run_sector_methods

mol = gto.M(
    atom="""
    F   0.0000   0.0000   0.0000
    H   0.0000   0.0000   0.9168
    """,
    unit='Angstrom', basis='6-31g', verbose=0,
)
mf = scf.RHF(mol).run(conv_tol=1e-12)

for sector in ('ip', 'ea'):
    results = run_sector_methods(
        mf, ['nD-ADC(3)', 'NRL3-ISR(3)'], frozen=1, sector=sector, nroots=3,
    )
    for method, poles in results.items():
        print(f'\n{method} {sector.upper()}')
        for root, pole in enumerate(poles):
            print(f'root {root}: {pole.binding_energy_ev:12.6f} eV, '
                  f'PS = {pole.strength:.6f}, residual = {pole.residual:.2e}')
            # cubegen.orbital(mol, 'dyson.cube', pole.dyson_ao) works unchanged.
