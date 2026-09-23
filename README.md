# Nondiagonal PyEPT: Reference Implementation of Non-diagonal Electron Propagator Methods

[License: MIT](LICENSE) · [Citation](CITATION.md) · [Contributing](CONTRIBUTING.md)

Non-diagonal electron propagator calculations in Python, using NumPy, SciPy, and PySCF. The original Dyson approximations follow the Opoku–Pawłowski–Ortiz articles cited below. Separate-sector ADC(3) and an explicitly experimental NRL3-derived representation are also available.

- E. Opoku, F. Pawłowski, J. V. Ortiz, *J. Chem. Phys.* **159**, 124109 (2023), [doi:10.1063/5.0168779](https://doi.org/10.1063/5.0168779).
- E. Opoku, F. Pawłowski, J. V. Ortiz, *J. Chem. Phys.* **155**, 204107 (2021), [doi:10.1063/5.0070849](https://doi.org/10.1063/5.0070849).

The core uses a **matrix-free Hermitian Hamiltonian**, explicit antisymmetric spin-orbital tensors, PySCF integrals, and residual-controlled Davidson iterations. The original Dyson methods retain both 2h1p and 2p1h manifolds. This is a research implementation with the validation scope described below.

## Implemented methods

| Name | Meaning |
|---|---|
| `ND2` / `ADC(2)` | Non-diagonal second-order Dyson self-energy |
| `2ph-TDA` | Two-particle-one-hole Tamm–Dancoff approximation |
| `NR2` | Non-diagonal renormalized second order |
| `NRP3` | Non-diagonal renormalized partial third order |
| `NRQ3` | Non-diagonal renormalized quasiparticle third order |
| `nD-NRL3` | Static opposite-sector extension of NRL3; [definition and validation](docs/NON_DYSON_NRL3.md) |
| `NRL3` | Non-diagonal renormalized linear third order |
| `ADC(3)` | Conventional Dyson ADC(3) with Schirmer–Angonoa DEM static self-energy |
| `3+`| Third-order plus |
| `BD-T1` | Brueckner-doubles reference with terms linear in doubles and triple operators |

`sector="ea"` applies the particle–hole counterpart of the asymmetric NR2/NRP3/NRQ3 truncations. ND2, 2ph-TDA, NRL3, 3+, and BD-T1 have the same Hamiltonian for IP and EA. Merely changing the sign of an IP does not implement an EA-specific NRQ3 calculation.

**Reference scope:** real, molecular, closed-shell RHF; BD-T1 automatically constructs a semicanonical Brueckner reference with PySCF BCCD. UHF, ROHF, DFT, density-fitted SCF references, complex/spinor orbitals, periodic systems and gradients are not implemented. Established non-Dyson ADC(3) is available through `SectorEPT` (PySCF >= 2.14); see below. The static `nD-NRL3` extension is supported; it is distinct from non-Dyson ADC. `ADC(3)` selects conventional Dyson ADC(3) with the DEM static correction; `ADC(3)-strict` continues to select `3+`. Full ADC(4), other static variants, and the articles' diagonal-only methods are outside this implementation.

## Install and run

Use Python 3.10 or newer in a virtual environment:

```bash
git clone https://github.com/ernopoku/Nondiagonal-PyEPT.git
cd Nondiagonal-PyEPT
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python examples/water.py
python -m nondiagonal_ept examples/hf.json -o hf_results.json
pytest -q
```

PySCF must have a wheel or a supported build on the chosen platform. Tested package versions are recorded in `docs/environment.txt`.

```python
from pyscf import gto, scf
from nondiagonal_ept import EPT

mol = gto.M(atom="H 0 0 0; F 0 0 0.9178", basis="cc-pvtz", verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)

calculation = EPT(mf, "NRL3", frozen=1, sector="ip")
for pole in calculation.kernel(targets=[4, 2], tol=1e-9):
    print(f"IP = {pole.binding_energy_ev:.6f} eV, PS = {pole.strength:.6f}, "
          f"residual = {pole.residual:.2e}")
    # pole.dyson_mo: active spatial-MO coefficients, norm squared = strength
    # pole.dyson_ao: AO coefficients, C.T @ S_AO @ C = strength
```

Targets are **zero-based indices in the ORIGINAL spatial MO list**, before freezing. For example, the fifth original spatial MO is always target `4`, including when a core orbital is frozen. `frozen=1` freezes the lowest occupied spatial MO; an explicit list can exclude occupied and/or virtual MOs. `spin=0` selects the alpha creator sector; `spin=1` selects beta. They have equal spectra for a closed-shell reference.

For attachment:

```python
attachment = EPT(mf, "NRQ3", frozen=1, sector="ea")
root = attachment.kernel(targets=[5])[0]
print(root.energy, root.binding_energy_ev)
```

`energy` is the signed propagator pole ω in hartree. For removal ω = E(N) − E(N−1), so IP = −ω. For addition ω = E(N+1) − E(N), so EA = −ω. Consequently, an unbound attachment has positive ω and negative EA. The program does not classify states using the sign of ω alone.

For BD-T1, use `EPT(mf, "BD-T1", frozen=1)`. This runs CCSD and orbital iterations, verifies the Brueckner singles norm, then builds the propagator. Target MO indices label the resulting semicanonical Brueckner orbitals; their character/order need not match the original canonical RHF orbitals. The supplied RHF object is preserved. Frozen core orbitals stay fixed at their input RHF coefficients throughout orbital optimization. For N2 use `frozen=2`; see `examples/n2_bdt1.py`.

## Solvers and observables

- `kernel(targets=[...])` follows the primary pole with largest overlap on each target simple operator. It returns only residual-converged roots or raises `ConvergenceError`. Strongly mixed primary/satellite states can require inspection of the full spectrum. Duplicate convergence from two targets raises an error.
- `dense_spectrum()` returns every pole for small systems, with a default dimension guard of 2500. It is also an independent eigensolver check on Davidson.
- `self_energy(E)` evaluates the non-diagonal Schur complement with checked MINRES solves. `self_energy(E, derivative=True)` returns its analytic energy derivative. Energies at triple-space singularities can fail to converge and raise an exception.
- A Dyson orbital is **not** unit-normalized by default. Its squared norm is its pole strength. `normalized_dyson_mo` provides the unit-normalized shape.
- The program uses spatial four-index ERIs and spin-block contractions, avoiding a full spin-orbital V⁴ tensor. It does not use density fitting or distributed/out-of-core propagator tensors. `max_memory_mb` is a conservative preflight estimate, not an operating-system memory limit. Very large basis sets remain expensive.

## References

- E. Opoku, F. Pawłowski, J. V. Ortiz, *J. Chem. Phys.* **159**, 124109 (2023), [doi:10.1063/5.0168779](https://doi.org/10.1063/5.0168779), especially Eqs. 18–24, 28–33 and Table IV.
- E. Opoku, F. Pawłowski, J. V. Ortiz, *J. Chem. Phys.* **155**, 204107 (2021), [doi:10.1063/5.0070849](https://doi.org/10.1063/5.0070849).
- An erratum exists at [doi:10.1063/5.0167154](https://doi.org/10.1063/5.0167154). 
- [PySCF AO-to-MO documentation](https://pyscf.org/contributor/ao2mo_developer.html). The BD reference uses `pyscf.cc.bccd.bccd_kernel_` and the spin-amplitude conversion in `pyscf.cc.addons`.

## Cartesian geometry and multiple methods

Use one Cartesian row per atom: element, x, y, z. Set `unit="Angstrom"`
explicitly. The following example uses F at the origin and a 0.9168 Å bond.
This differs from the earlier 0.9178 Å example, so energies will differ slightly.

```python
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
              f"PS = {pole.strength:.6f}, residual = {pole.residual:.2e}")
```

`run_methods` reuses the converged RHF reference and runs each requested method
sequentially. It returns a dictionary keyed by canonical method name, with a
list of poles for each method. Propagator intermediates are built separately
for each method. All EPT options and solver settings are shared. Duplicate
methods (including aliases) are rejected; convergence failures raise an exception.
The existing single-method `EPT` interface remains available.


## License

Nondiagonal PyEPT is distributed under the [MIT license](LICENSE).
Copyright (c) 2026 Ernest Opoku. Dependencies retain their own licenses.

## Citation

Please cite the software and the relevant scientific method papers when using
this work in research. See [citation guidance and BibTeX](CITATION.md), or use
GitHub's **Cite this repository** control, powered by [CITATION.cff](CITATION.cff).
Record the exact Git commit and computational settings for reproducibility.

## Contributions

Bug reports, method improvements, tests, examples, and documentation contributions
are welcome. Read [CONTRIBUTING.md](CONTRIBUTING.md) for setup, numerical validation,
and pull-request guidance. Report problems through
[GitHub Issues](https://github.com/ernopoku/Nondiagonal-PyEPT/issues).

Printed examples label pole strength as **PS**. The Python attribute `pole.strength` and machine-readable JSON fields retain their existing names for compatibility.

## Static non-Dyson NRL3

Use `EPT(mf, "nD-NRL3", frozen=1, sector="ip")` or include `"nD-NRL3"` in `run_methods`. For IP, 2ph contributions are frozen at HF orbital energies and symmetrized; for EA, 2hp contributions are frozen instead. See [equations, PS interpretation, examples, and validation limits](docs/NON_DYSON_NRL3.md).


## Separate IP/EA intermediate-state representations

```python
from nondiagonal_ept import SectorEPT, run_sector_methods

# Established non-Dyson ADC(3); counts states rather than targeting MO indices.
poles = SectorEPT(mf, "nD-ADC(3)", frozen=1, sector="ip").kernel(nroots=3)

# For a SMALL basis: compare with the experimental NRL3-derived sector model.
results = run_sector_methods(mf, ["nD-ADC(3)", "NRL3-ISR(3)"],
                             frozen=1, sector="ip", nroots=3)
```

Install `python -m pip install -e '.[sector]'` for the ADC dependency.
`NRL3-ISR(3)` is a third-order canonical reduction of the NRL3 auxiliary matrix,
with consistently transformed Dyson amplitudes. It is **experimental**, uses
dense matrices (default full-dimension limit 1200), and is not a separately
derived ground-state-metric non-Dyson NRL3 theory. It is distinct from the
static `nD-NRL3` extension. No large-basis speedup is claimed.

Read the [derivation, API, PS conventions and validation](docs/SECTOR_ISR.md)
or the [step-by-step Wiki guide](https://github.com/ernopoku/Nondiagonal-PyEPT/wiki/Separate-Sector-Methods).
Run `python examples/sector_methods.py` for a small working example.


## Conventional Dyson ADC

```python
results = run_methods(mf, ["ADC(2)", "ADC(3)", "3+"],
                      frozen=1, sector="ip", targets=[4, 2])
```

`ADC(2)` is the existing `ND2` method (its result key remains `ND2`).
`ADC(3)` now includes the Schirmer–Angonoa Dyson-expansion static self-energy,
correct through fourth order with selected higher-order contributions. The
dynamic blocks retain third-order ADC accuracy. This differs from both the
strict static `3+` method and `SectorEPT(mf, "nD-ADC(3)")`.

Use `static_tol=1e-10` and `static_max_cycle=500` to control the static solves;
`kernel` controls remain separate. All solves are residual checked. See
[working equations, examples and validation](docs/DYSON_ADC.md),
[the Wiki guide](https://github.com/ernopoku/Nondiagonal-PyEPT/wiki/Dyson-ADC),
and `examples/dyson_adc.py`. The JSON CLI also accepts these methods.
