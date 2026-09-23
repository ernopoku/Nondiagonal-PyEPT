# Conventional Dyson ADC(2) and ADC(3)

Use `EPT` for these methods. Both IP and EA are contained in the same Dyson Hamiltonian.

| Name | Definition |
|---|---|
| `ADC(2)` / `ND2` | Conventional second-order Dyson self-energy |
| `ADC(3)` / `ADC(3)-DEM` | Conventional Dyson ADC(3) with Schirmer–Angonoa DEM static self-energy |
| `ADC(3)-strict` / `3+` | Strict third-order static self-energy, with the same ADC(3) dynamic blocks |
| `nD-ADC(3)` | Separate IP/EA non-Dyson ADC(3), available through `SectorEPT` |

`ADC(2)` was already implemented as ND2; its canonical result key remains `ND2`. `ADC(3)` now selects the DEM formulation rather than raising an ambiguity error. Existing `3+` and non-Dyson results are unchanged.

## Update

```bash
git pull
python -m pip install -e .
```

## Calculate several methods

```python
from pyscf import gto, scf
from nondiagonal_ept import EPT, run_methods

mol = gto.M(
    atom="""
    F   0.0000   0.0000   0.0000
    H   0.0000   0.0000   0.9168
    """,
    unit="Angstrom", basis="cc-pvtz", verbose=0,
)
mf = scf.RHF(mol).run(conv_tol=1e-12)

results = run_methods(
    mf, ["ADC(2)", "ADC(3)", "3+"],
    frozen=1, sector="ip", targets=[4, 2], tol=1e-9,
)
for method, poles in results.items():
    for pole in poles:
        print(f"{method} MO {pole.target}: IP = {pole.binding_energy_ev:.6f} eV, "
              f"PS = {pole.strength:.6f}, residual = {pole.residual:.2e}")
```

For attachment:

```python
calc = EPT(mf, "ADC(3)", frozen=1, sector="ea")
for pole in calc.kernel(targets=[5]):
    print(f"EA = {pole.binding_energy_ev:.6f} eV, PS = {pole.strength:.6f}")
```

Targets are zero-based original spatial MO indices, not root ranks. `pole.energy` is signed omega in hartree; binding IP/EA is minus omega. An unbound attachment therefore has negative binding EA. The existing per-spin PS, MO/AO Dyson orbitals and [cube visualization](Visualizing-Dyson-Orbitals) are supported.

## Static convergence controls

```python
calc = EPT(mf, "ADC(3)", frozen=1,
           static_tol=1e-10, static_max_cycle=500)
print(calc.hamiltonian.static_diagnostics["static_residual"])
poles = calc.kernel(targets=[4, 2], tol=1e-9)
```

The static correction is computed once during construction. Triple-space resolvents and the occupied–virtual response equation use matrix-free MINRES, with explicit residual checks. `static_max_cycle` applies per static solve; the eigenvalue solver has its own `max_cycle/max_space`. A failed solve raises `ConvergenceError`, without silently substituting `3+`.

The JSON CLI also accepts both methods and these controls:

```bash
python -m nondiagonal_ept examples/hf_adc3.json -o hf_adc3_results.json
python examples/dyson_adc.py
```

## Meaning of the ADC(3) static correction

The dynamic vertices and triple blocks are the same as in `3+`. The DEM static matrix instead solves a linear equation for the correlated-density correction, holding the RHF orbitals and dynamic blocks fixed. The constant self-energy is correct through fourth order with selected higher-order terms; this does not make the full method ADC(4). It is distinct from a fully self-consistent Dyson-density iteration and from non-Dyson Sigma(4+).

See [the complete equations and references](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/DYSON_ADC.md). The DEM correlation-density trace need not vanish exactly. No particle-number rescaling is applied; diagnostics expose its spin-summed trace.

## Validation

All **83 tests pass**. Independent checks include contour integration, a full spin-orbital linear system, fourth-order static accuracy against small-system FCI densities, IP/EA residues, frozen virtual orbitals, phase/spin invariance and failed-solver handling.

39 molecular poles across HF/cc-pVDZ, HF/cc-pVTZ, water/cc-pVDZ and N2/cc-pVDZ have residuals below 1e-9 hartree. These are internal comparisons, not validation against an independent molecular ADC(3)-DEM program. [Raw data](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/dyson_adc_comparison.json) can be regenerated with `python examples/validate_dyson_adc.py`.

Real molecular canonical closed-shell RHF is required. Frozen orbitals retain their mean-field contribution; the correlation response uses active orbitals only. Large calculations remain memory intensive, although the static solver does not build a dense full propagator matrix.

[Home](Home) · [Separate non-Dyson methods](Separate-Sector-Methods) · [Citation guidance](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/CITATION.md)
