# Validation report

Validation performed on 22 September 2026 using Python 3.12, NumPy, SciPy, and PySCF 2.14.0 on macOS arm64. Exact installed versions are recorded in `environment.txt`. All reported pole energies are in hartree unless stated otherwise.

## What was checked

The automated tests cover:

1. First- and second-order simple–triple vertices against explicit fermionic operator matrices.
2. Both first-order triple blocks against independent N−1/N+1 determinant Hamiltonians.
3. Linear doubles corrections to both triple blocks against the Hermitized super-operator metric.
4. Static third-order densities and second-order singles against reference-state and Rayleigh–Schrödinger calculations.
5. Strict 3+ accuracy through third order: halving the fluctuation potential reduces the leading full-CI error by approximately 16, as required for a fourth-order remainder.
6. PySCF MP2 correlation-energy agreement with the spin-orbital amplitude convention.
7. Davidson roots against full dense diagonalization for all seven HF-reference methods.
8. Pole-strength sum rules, AO and MO Dyson normalization, alpha/beta equality, and orbital-phase invariance.
9. Analytic self-energy derivatives and the residue derivative identity.
10. Invalid references, frozen targets, ambiguous ADC names, and deliberately unconverged solver calls.
11. Brueckner preparation, preservation of the caller's reference, matrix Hermiticity, and diagonal preconditioning.

The final pytest output is saved in `test_results.txt`. Determinant-algebra tests include both two- and four-electron references; the oracle does not use the production einsum equations.

## Comparison with supplied reference runs

The reproducible script `examples/validate_reference.py` uses the supplied geometries, spherical cc-pVTZ basis, and frozen-core counts. It compares NRP3, NRQ3, and NRL3 for IP and EA roots of HF, H2O, N2, F2, and CO. Both extra N2 NRL3 EA roots in the archive are included.

**58 poles were compared.** All Python eigenpair residuals are below 1.895e-09 hartree. 53 of 58 archived energies agree within 1e-5 hartree (0.000272 eV). The five larger differences are diagnosed below. Pole strengths are reported without filtering; they can be more sensitive to incomplete eigenvector convergence than energies.

Raw comparisons, including every energy, strength, residual, iteration count, and source filename, are in `reference_comparison.json`. `tests/reference_poles.json` contains all 76 extracted reference poles, including the ethylene cases not executed in this validation run. The reference JSON is not a promise that the legacy stopping criterion was sufficient.

The PySCF and archived SCF energies differ slightly (approximately 1e-7 hartree in these comparisons). Basis-library numerical differences and SCF thresholds can therefore also affect sub-microhartree comparisons. We do not require bitwise equality between different integral engines.

## Confirmed legacy premature convergence

The archived reference calculations used successive Ritz energy changes after at least three iterations as their stopping criterion, without an eigenvector-residual check.

To distinguish an equation error from an iteration error, we reproduced the legacy H0-preconditioned Davidson iterations on the Python Hamiltonian, using overlap tracking. In each of the five discrepant cases, the legacy stopping iteration reproduces the archived energy within about 1.3e-7 hartree and recovers its pole strength, but leaves a large residual. Continuing the same equations removes the false convergence.

| Molecule / method | Original MO (0-based) | Legacy stopping iteration | reference pole | Reproduced early pole | Residual at that iteration |
|---|---:|---:|---:|---:|---:|
| H2O / NRL3 | 3 | 3 | -0.5381424623 | -0.5381425596 | 0.171611 |
| N2 / NRL3 | 3 | 4 | -0.6858683991 | -0.6858685190 | 0.198708 |
| CO / NRP3 | 6 | 5 | -0.5161996810 | -0.5161997397 | 0.158753 |
| CO / NRP3 | 5 | 3 | -0.6210800406 | -0.6210801596 | 0.221665 |
| CO / NRQ3 | 5 | 3 | -0.6240965037 | -0.6240966157 | 0.229715 |

The production code deliberately does **not** reproduce those unconverged answers.

| Molecule / method | MO | Archived pole | Residual-converged pole | Pole change (eV) |
|---|---:|---:|---:|---:|
| H2O / NRL3 | 3 | -0.5381424623 | -0.5374795285 | +0.018039 |
| N2 / NRL3 | 3 | -0.6858683991 | -0.6836542530 | +0.060250 |
| CO / NRP3 | 6 | -0.5161996810 | -0.5162355768 | -0.000977 |
| CO / NRP3 | 5 | -0.6210800406 | -0.6186947408 | +0.064907 |
| CO / NRQ3 | 5 | -0.6240965037 | -0.6216053119 | +0.067789 |

The last column is the change in signed propagator pole; the corresponding IP change has the opposite sign. The largest change is approximately 0.068 eV. The comparison is evidence of a solver convergence defect in these archived runs, not evidence that the physical accuracy of the approximation improved by the same amount.

Each trace is saved beside this report. Reproduce examples with:

```bash
python examples/diagnose_convergence.py --molecule H2O --method NRL3 --target 3
python examples/diagnose_convergence.py --molecule N2 --method NRL3 --target 3
python examples/diagnose_convergence.py --molecule CO --method NRQ3 --target 5
```

## BD-T1 validation boundary

BD-T1 uses converged PySCF Brueckner doubles, half-weight linear vertex corrections, the off-diagonal Brueckner Fock block, and the linear doubles triple-block corrections. Its tensor equations are checked against determinant operator algebra. H2/6-31G and frozen-core H2O/STO-3G execute successfully; the latter's values and Brueckner singles norm are in `bdt1_smoke.json`.

The supplied N2/cc-pVDZ fixed-core reference now provides a molecular regression for three BD-T1 IP poles. See [BD-T1 validation](BD_T1_VALIDATION.md) for the correction, comparisons, and precision limits. ND2, 2ph-TDA, NR2, and strict 3+ retain their independent algebra and molecular checks.

## Remaining scope limits

- The 2021 erratum was identified but its full text was unavailable. New-method block definitions were taken from the supplied 2023 article and checked against operator algebra and archived numerical results.
- This validates the stated restricted-reference methods, not UHF/ROHF, fourth-order static ADC variants, all diagonal methods, experimental order-flag combinations, or the articles' full benchmark datasets.
- The current implementation is an in-core research code with matrix-free propagator iterations. Performance has not been characterized for very large molecules or distributed execution.

## Reproduction

From the package directory, after installation:

```bash
pytest -q
python examples/validate_reference.py --molecules HF H2O N2 F2 CO --output comparison.json
python -m nondiagonal_ept examples/hf.json -o hf_results.json
```

Some BLAS builds run small tensor contractions faster with one thread. The recorded benchmark used `OPENBLAS_NUM_THREADS=1` and `VECLIB_MAXIMUM_THREADS=1`; scientific results must not depend on those settings beyond rounding.

## Static nD-NRL3 extension

See [the dedicated validation and definition](NON_DYSON_NRL3.md) and [ten IP/EA comparisons](non_dyson_comparison.json). This extension has independent dense-resolvent and weak-coupling checks but no external nD-NRL3 benchmark implementation. The full suite now passes 55 tests.
