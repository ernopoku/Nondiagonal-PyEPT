# Sector-projected nD-NRL3: cross-basis validation

These data test the explicitly revised `sector-projected-v2` construction. The [method guide](NON_DYSON_NRL3.md) gives the working equations, energy-order argument, and limitations. This comparison uses NRL3 as the parent method; it does not establish exact or experimental accuracy.

## Protocol

All cases use neutral singlet RHF, Angstrom coordinates, RHF tolerance 1e-12, frozen lowest 1s orbitals on non-hydrogen atoms, static residual tolerance 1e-10, and eigenpair residual tolerance 1e-9. No virtual orbitals are frozen. Original zero-based MO indices are recorded. The comparison solvers allow 300 iterations and a 60-vector subspace.

HF, H2O, N2, and CO are tested with cc-pVDZ, aug-cc-pVDZ, cc-pVTZ, and aug-cc-pVTZ. HF additionally uses cc-pVQZ and aug-cc-pVQZ. User-supplied BeO, Li2, LiF, LiH, and C3 geometries are tested with the four DZ/TZ bases. The supplied set includes the highest two active occupied targets where available and the lowest two virtual targets; LiH and Li2 have only one active occupied target after freezing.

RHF internal and RHF-to-UHF external stability checks are recorded for the supplied set. Internally unstable references are restarted at most three times; unresolved internal instability stops that case. External instability is reported while preserving the RHF reference, since these implementations do not support UHF. The core set did not undergo that additional stability screen.

During validation, residual checking was optimized to omit interactions in the unused triple sector on identically zero input. The independent full-parent matrix tests still pass; this changes runtime, not the operator. The record retains the pre-optimization source hash and per-case hashes when available, alongside the final source hashes. The pending large HF EA calculation was restarted after this optimization; no numerical failure was reported for the interrupted attempt.

The first C3/aug-cc-pVTZ attempt stopped at the 16000 MB workspace guard (estimated need 21302 MB). It was retried separately at 24000 MB on a machine with 32 GiB physical RAM; the record preserves both attempts. This is a resource limit, not a convergence failure. Other cases used 16000 MB.

Core calculations requested eight threads; the supplied set requested one. The two sweeps overlapped in time, so elapsed values are diagnostics, not controlled performance comparisons. The environment and source hashes are included in the [machine-readable record](non_dyson_sector_validation.json).

There are **130 completed pole comparisons** in this record. The full automated test suite passes **153 tests**, including an HF/cc-pVTZ public-API regression, legacy energy reproduction, independent dense matrix tests, both spins/sectors, fourth-order weak-coupling differences, energy-origin invariance, original-MO mapping, derivatives, PS/AO normalization, and residual checks. The JSON CLI was also exercised with explicit `static_space`. The molecular sweeps are additional to these tests.

## Maximum absolute differences from NRL3

Each entry is the maximum over the requested targets in that molecule/basis/sector. A missing completed sector is shown as a dash; see failures below. Values are in eV.

| Molecule | Basis | Max IP difference | Max EA difference |
|---|---|---:|---:|
| HF | cc-pvdz | 0.043719 | 0.004097 |
| HF | aug-cc-pvdz | 0.061175 | 0.000983 |
| HF | cc-pvtz | 0.047996 | 0.004494 |
| HF | aug-cc-pvtz | 0.054394 | 0.000929 |
| HF | cc-pvqz | 0.049302 | 0.004115 |
| HF | aug-cc-pvqz | 0.051644 | 0.000903 |
| H2O | cc-pvdz | 0.041855 | 0.004226 |
| H2O | aug-cc-pvdz | 0.059705 | 0.000980 |
| H2O | cc-pvtz | 0.045582 | 0.004432 |
| H2O | aug-cc-pvtz | 0.052076 | 0.000890 |
| N2 | cc-pvdz | 0.075787 | 0.027413 |
| N2 | aug-cc-pvdz | 0.083624 | 0.000511 |
| N2 | cc-pvtz | 0.077227 | 0.042733 |
| N2 | aug-cc-pvtz | 0.078864 | 0.000426 |
| CO | cc-pvdz | 0.049076 | 0.033463 |
| CO | aug-cc-pvdz | 0.051136 | 0.016680 |
| CO | cc-pvtz | 0.046569 | 0.041195 |
| CO | aug-cc-pvtz | 0.046216 | 0.001003 |
| LiH | cc-pvdz | 0.022517 | 0.000704 |
| LiH | aug-cc-pvdz | 0.022655 | 0.000403 |
| LiH | cc-pvtz | 0.019575 | 0.000548 |
| LiH | aug-cc-pvtz | 0.019846 | 0.000401 |
| Li2 | cc-pvdz | 0.009894 | 0.011183 |
| Li2 | aug-cc-pvdz | 0.010660 | 0.007343 |
| Li2 | cc-pvtz | 0.008097 | 0.010809 |
| Li2 | aug-cc-pvtz | 0.008562 | 0.006261 |
| LiF | cc-pvdz | 0.056899 | 0.001851 |
| LiF | aug-cc-pvdz | 0.065288 | 0.000162 |
| LiF | cc-pvtz | 0.057424 | 0.000805 |
| LiF | aug-cc-pvtz | 0.060276 | 0.000157 |
| BeO | cc-pvdz | 0.073235 | 0.022226 |
| BeO | aug-cc-pvdz | 0.071980 | 0.017734 |
| BeO | cc-pvtz | 0.069112 | 0.018228 |
| BeO | aug-cc-pvtz | 0.068390 | 0.016674 |
| C3 | cc-pvdz | 0.019349 | 0.046062 |
| C3 | aug-cc-pvdz | 0.016426 | 0.054380 |
| C3 | cc-pvtz | 0.014425 | 0.058271 |
| C3 | aug-cc-pvtz | 0.013639 | 0.061330 |

## Original HF discrepancy

| Basis | MO | Legacy nD − NRL3 (eV) | Version 2 − NRL3 (eV) |
|---|---:|---:|---:|
| cc-pvdz | 4 | -0.032836 | -0.043719 |
| cc-pvdz | 2 | -0.010217 | -0.032762 |
| aug-cc-pvdz | 4 | -0.038651 | -0.061175 |
| aug-cc-pvdz | 2 | +0.041347 | -0.040132 |
| cc-pvtz | 4 | +0.602018 | -0.047996 |
| cc-pvtz | 2 | +0.450276 | -0.034975 |
| aug-cc-pvtz | 4 | +0.270351 | -0.054394 |
| aug-cc-pvtz | 2 | +2.783286 | -0.036545 |

Version 2 removes the observed wrong-sector sampling failure rather than fitting these energy differences. No case-specific shifts or tolerance adjustments were introduced. Remaining differences are approximation differences, not expected to vanish exactly.

## Reference stability, root matching, and failures

- LiH: RHF-to-UHF external instability in none of the completed stability checks.
- Li2: RHF-to-UHF external instability in cc-pvdz, aug-cc-pvdz, cc-pvtz, aug-cc-pvtz.
- LiF: RHF-to-UHF external instability in none of the completed stability checks.
- BeO: RHF-to-UHF external instability in none of the completed stability checks.
- C3: RHF-to-UHF external instability in cc-pvdz, aug-cc-pvdz, cc-pvtz, aug-cc-pvtz.

External instability limits the physical interpretation of these RHF-based results, even if nD-NRL3 and NRL3 agree. It is not fixed by projecting the propagator space.

All completed paired roots have normalized Dyson overlap squared above 0.9.

All final case attempts completed without a convergence failure; the resolved C3 memory-guard stop is recorded above.

## Reproduce

```sh
python -m pytest -q
python examples/validate_non_dyson_core.py
python examples/validate_non_dyson_sector.py --threads 1 --memory-mb 24000 --output supplied-validation
```

The supplied-set runner saves each stage, per-root results, and errors, and uses a 900-second per-case timeout by default (`--timeout` changes it). The separate C3 retry allowed 1200 seconds. Failed/incomplete cases produce a nonzero exit status in the released runner. Both scripts record actual geometries and frozen spaces. Full raw per-case static diagnostics are generated by the scripts; the repository record retains maxima instead of every iteration/residual entry.

For legacy reproduction set `static_space="full"`; the historical audit script does this explicitly. The ten-case comparison file `non_dyson_comparison.json` is regenerated for version 2, while `non_dyson_comparison_v1.json` preserves the original results.

No finite molecular benchmark guarantees accuracy across all basis sets or strongly correlated systems. In particular, opposite-sector resonances within the retained sampling range, weak gaps, reference instability, or strong state mixing remain reasons to inspect results carefully. Old full-space manuscript formulas, tables, dimensions, and amplitudes must not be presented as version-2 results.
