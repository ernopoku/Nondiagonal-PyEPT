# Separate IP/EA methods

Two methods are available through `SectorEPT`:

| Name | What it computes | Status |
|---|---|---|
| `nD-ADC(3)` | Established non-Dyson IP/EA ADC(3), using PySCF | Checked against direct PySCF |
| `NRL3-ISR(3)` | Third-order canonical sector reduction of the NRL3 auxiliary matrix | Experimental, dense small-system prototype |

**NRL3-ISR(3) is not a separately derived ground-state-metric non-Dyson NRL3 theory.** It has explicit NRL3-specific matrix and transition-moment equations, but no independent published benchmark. It is distinct from both ADC(3) and the [static nD-NRL3 extension](Non-Dyson-NRL3).

## 1. Update and install

From your repository directory:

```bash
git pull
python -m pip install -e ".[sector]"
```

The ADC adapter requires PySCF 2.14 or newer. Both interfaces currently support real, canonical, molecular closed-shell RHF and conventional integrals.

## 2. Run either method, or both

```python
from pyscf import gto, scf
from nondiagonal_ept import SectorEPT, run_sector_methods

mol = gto.M(
    atom="""
    F   0.0000   0.0000   0.0000
    H   0.0000   0.0000   0.9168
    """,
    unit="Angstrom", basis="6-31g", verbose=0,
)
mf = scf.RHF(mol).run(conv_tol=1e-12)

results = run_sector_methods(
    mf, ["nD-ADC(3)", "NRL3-ISR(3)"],
    frozen=1, sector="ip", nroots=3, tol=1e-9,
)
for method, poles in results.items():
    for root, pole in enumerate(poles):
        print(f"{method} root {root}: IP = {pole.binding_energy_ev:.6f} eV, "
              f"PS = {pole.strength:.6f}, residual = {pole.residual:.2e}")
```

For a single method, including attachment:

```python
calc = SectorEPT(mf, "nD-ADC(3)", frozen=1, sector="ea")
poles = calc.kernel(nroots=3)
for root, pole in enumerate(poles):
    print(f"root {root}: EA = {pole.binding_energy_ev:.6f} eV, PS = {pole.strength:.6f}")
```

Use a larger basis such as cc-pVTZ with `nD-ADC(3)` alone when your resources permit. The small basis above allows the experimental dense method to run too; it is not a basis recommendation for accurate affinities.

## 3. Understand root labels and signs

`nroots` counts states, not MO indices. IP roots are ordered by increasing ionization cost; EA roots by increasing attachment energy (largest binding EA first). Low-PS satellites can appear. Compare Dyson character when matching states between methods. Degenerate Dyson orbitals may rotate within their degenerate subspace.

`pole.target` is `None`. This API uses `SectorEPT.kernel(nroots=...)`, not `EPT.kernel(targets=...)` or the existing JSON CLI. `frozen` retains the original count/index-list convention. `calc.original_mos` maps active Dyson MO coefficients to original spatial MO indices.

`pole.energy` is signed omega in hartree; `binding_energy_ev` is minus omega in eV. A bound IP is positive; an unbound attachment has negative binding EA. PS is reported for one spin. The restricted PySCF spectroscopic factor includes both spins, so this adapter returns half that factor. Its Dyson amplitudes already have the correct one-spin normalization.

## 4. Visualize the Dyson orbital

```python
from pyscf.tools import cubegen
pole = results["nD-ADC(3)"][0]
cubegen.orbital(mol, "ip_root0.cube", pole.dyson_ao,
                nx=100, ny=100, nz=100, margin=5.0)
```

Open the cube in your molecular viewer and display positive and negative isosurfaces. The MO squared norm and AO overlap-metric norm both equal PS. These orbitals use the full transition matrix, not just the simple part of the eigenvector. Follow [Visualizing Dyson Orbitals](Visualizing-Dyson-Orbitals) for detailed viewing instructions.

## 5. Experimental NRL3 derivation and limits

NRL3-ISR(3) partitions occupied simple plus 2hp configurations into IP and virtual simple plus 2ph into EA. An antisymmetric canonical generator removes cross-sector coupling through third order. Separate Hermitian sector matrices are then diagonalized. The same transformation supplies transition amplitudes; an explicitly documented exponential completion preserves positivity and combined spectral completeness.

Read the [full derivation and transition-moment equations](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/SECTOR_ISR.md). A direct many-electron derivation from the NRL3 Hermitized metric remains unestablished. Size consistency and strong-correlation behavior of this new completion have not been established.

The prototype builds the full auxiliary matrix before reducing it. It offers no current speed or memory advantage over matrix-free NRL3. Its default full-dimension limit is 1200; larger cc-pVTZ examples usually exceed this. Use small systems first. `max_cycle` and `max_space` affect ADC only.

```python
calc = SectorEPT(mf, "NRL3-ISR(3)", frozen=1, sector="ip")
print(calc.diagnostics)
```

Diagnostics report denominator gaps, generator norms and omitted higher-order/cross-sector terms. They are not rigorous error bounds. Eigenpair residuals verify the numerical solve, not physical accuracy.

## 6. Validation and reproduction

All 70 tests pass. ADC is compared against direct PySCF in both sectors, including frozen virtual orbitals, transition amplitudes and per-spin PS. The NRL3 construction has independent matrix-exponential and fourth-order remainder checks for energies, propagators and Dyson residues, plus zero-interaction and spectral-completeness tests.

In 18 small-basis IP/EA comparisons for HF, water and N2, NRL3-ISR(3) differs from parent NRL3 by at most 0.024 eV. This is an internal consistency comparison, not independent validation of a new non-Dyson theory. Some included roots are nearly dark satellites.

```bash
python examples/sector_methods.py
python examples/validate_sector_methods.py
pytest -q
```

[Raw comparisons](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/sector_comparison.json) · [Citations](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/CITATION.md) · [Home](Home)
