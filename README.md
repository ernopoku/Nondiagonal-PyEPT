# Nondiagonal PyEPT

[License: MIT](LICENSE) · [Citation](CITATION.md) · [Contributing](CONTRIBUTING.md)

Non-diagonal electron propagator calculations in Python, using NumPy, SciPy, and PySCF. The implemented approximations follow the Opoku–Pawłowski–Ortiz articles cited below.

The core uses a **matrix-free Hermitian Hamiltonian**, explicit antisymmetric spin-orbital tensors, PySCF integrals, and residual-controlled Davidson iterations. Both 2h1p and 2p1h manifolds are retained. This is a research implementation with the validation scope described below.

## Implemented methods

| Name | Meaning |
|---|---|
| `ND2` / `ADC(2)` | Non-diagonal second-order Dyson self-energy |
| `2ph-TDA` | First-order interactions in both triple manifolds |
| `NR2` | Non-diagonal renormalized second order |
| `NRP3` | Non-diagonal renormalized partial third order |
| `NRQ3` | Non-diagonal renormalized quasiparticle third order |
| `NRL3` | Non-diagonal renormalized linear third order |
| `3+` / `ADC(3)-strict` | Strict third-order Dyson ADC with ring/ladder renormalization |
| `BD-T1` | Brueckner-doubles reference with terms linear in doubles and triple operators |

`sector="ea"` applies the particle–hole counterpart of the asymmetric NR2/NRP3/NRQ3 truncations. ND2, 2ph-TDA, NRL3, 3+, and BD-T1 have the same Hamiltonian for IP and EA. Merely changing the sign of an IP does not implement an EA-specific NRQ3 calculation.

**Reference scope:** real, molecular, closed-shell RHF; BD-T1 automatically constructs a semicanonical Brueckner reference with PySCF BCCD. UHF, ROHF, DFT, density-fitted SCF references, complex/spinor orbitals, periodic systems, gradients, and non-Dyson ADC are not implemented. The generic name `ADC(3)` is deliberately rejected: use `ADC(3)-strict`. Fourth-order/renormalized-static ADC(3) variants and the articles' diagonal-only methods are outside this implementation.

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
    print(pole.binding_energy_ev, pole.strength, pole.residual)
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

For BD-T1, use `EPT(mf, "BD-T1", frozen=1)`. This runs CCSD and orbital iterations, verifies the Brueckner singles norm, then builds the propagator. Target MO indices label the resulting semicanonical Brueckner orbitals; their character/order need not match the original canonical RHF orbitals. The supplied RHF object is preserved.

## Solvers and observables

- `kernel(targets=[...])` follows the primary pole with largest overlap on each target simple operator. It returns only residual-converged roots or raises `ConvergenceError`. Strongly mixed primary/satellite states can require inspection of the full spectrum. Duplicate convergence from two targets raises an error.
- `dense_spectrum()` returns every pole for small systems, with a default dimension guard of 2500. It is also an independent eigensolver check on Davidson.
- `self_energy(E)` evaluates the non-diagonal Schur complement with checked MINRES solves. `self_energy(E, derivative=True)` returns its analytic energy derivative. Energies at triple-space singularities can fail to converge and raise an exception.
- A Dyson orbital is **not** unit-normalized by default. Its squared norm is its pole strength. `normalized_dyson_mo` provides the unit-normalized shape.
- The program uses spatial four-index ERIs and spin-block contractions, avoiding a full spin-orbital V⁴ tensor. It does not use density fitting or distributed/out-of-core propagator tensors. `max_memory_mb` is a conservative preflight estimate, not an operating-system memory limit. Very large basis sets remain expensive.

## Validation and a legacy convergence defect

See [the validation report](docs/VALIDATION.md), [equations](docs/THEORY.md), and the [project Wiki](https://github.com/ernopoku/Nondiagonal-PyEPT/wiki).

Independent tests evaluate fermionic operators in determinant spaces, compare MP2 amplitudes with PySCF, check the third-order limit against full CI, compare Davidson with dense diagonalization, and verify spectral sum rules, spin symmetry, MO-phase invariance, Dyson normalization, and analytic self-energy derivatives.

The archived reference results are valuable references but are **not all converged eigenpairs**. In particular, their water NRL3 calculation stops after three iterations because the energy change is small although the eigenvector residual is large. Reproducing those iterations recovers the archived energy and pole strength; continuing them changes the result. The rewrite fixes this by requiring a small residual. The report distinguishes this intentional numerical correction from agreement with the archived output.

BD-T1 is implemented and algebraically tested, with small-molecule execution tests. No independent BD-T1 molecular benchmark was supplied, so its quantitative validation is less extensive than NRP3/NRQ3/NRL3. No claim is made to reproduce the complete articles' statistical benchmark datasets.

## References

- E. Opoku, F. Pawłowski, J. V. Ortiz, *J. Chem. Phys.* **159**, 124109 (2023), [doi:10.1063/5.0168779](https://doi.org/10.1063/5.0168779), especially Eqs. 18–24, 28–33 and Table IV.
- E. Opoku, F. Pawłowski, J. V. Ortiz, *J. Chem. Phys.* **155**, 204107 (2021), [doi:10.1063/5.0070849](https://doi.org/10.1063/5.0070849).
- An erratum exists at [doi:10.1063/5.0167154](https://doi.org/10.1063/5.0167154). Its full text was not available during this work; the implemented new-method definitions use the supplied 2023 article and explicitly documented block definitions, rather than assuming the uncorrected 2021 formulas are definitive.
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
              f"strength = {pole.strength:.6f}, residual = {pole.residual:.2e}")
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
