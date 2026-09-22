from pyscf import gto, scf
from nondiagonal_ept import run_methods

mol = gto.M(
    atom="""
    F   0.0000   0.0000   0.0000
    H   0.0000   0.0000   0.9168
    """,
    unit="Angstrom",
    basis="cc-pvtz",
    verbose=0,
)
mf = scf.RHF(mol).run(conv_tol=1e-12)

results = run_methods(
    mf, ["NRL3", "NRQ3", "NRP3"],
    frozen=1, sector="ip", targets=[4, 2], tol=1e-9,
)
for method, poles in results.items():
    print(f"\n{method}")
    for pole in poles:
        print(f"MO {pole.target}: IP = {pole.binding_energy_ev:.6f} eV, "
              f"strength = {pole.strength:.6f}, residual = {pole.residual:.2e}")
