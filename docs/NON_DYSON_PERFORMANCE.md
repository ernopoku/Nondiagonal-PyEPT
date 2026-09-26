# nD-NRL3 performance and numerical equivalence

**Historical report:** these measurements concern `full-space-v1`. The default is now the sector-projected version 2; see the [current method guide](NON_DYSON_NRL3.md) and [new validation](NON_DYSON_SECTOR_VALIDATION.md).

The static opposite-sector construction is unchanged. This update accelerates its linear solves and contractions; it does not substitute a diagonal opposite-sector resolvent or remove virtual simple orbitals.

## Why the original implementation was slow

NRL3 iterates a few selected poles. nD-NRL3 first solves a shifted opposite-sector problem for **every active simple orbital**, including virtual orbitals in an IP calculation. The original implementation used unpreconditioned MINRES, evaluated the unused sector on every product, and repeatedly repacked quartic virtual integrals. Tiny residual-refinement problems could also be oversolved by many orders of magnitude. High virtual shifts may lie within the opposite-sector spectrum and are especially demanding.

## Changes

- Positive-definite, clipped diagonal preconditioning, including indefinite shifted systems.
- Exact products confined to the eliminated auxiliary sector.
- Packed nonzero spin blocks, with antisymmetry used to reconstruct the mixed-spin partner.
- Symmetric/antisymmetric virtual-pair kernels in reusable BLAS layout, or streamed slabs when the cache budget is insufficient.
- Periodic checks of the actual residual stop unnecessary iterations, particularly during refinement.
- Explicit static tolerance and iteration controls, progress messages, and JSON-serializable diagnostics.

The full shifted-operator residual threshold remains `1e-10 * max(1, norm(rhs))` by default. The preconditioner floor never changes the physical denominator. See [the method guide](NON_DYSON_NRL3.md#static-setup-controls-and-progress) for formulas, controls, and memory details.

## Measured timings

Runs used the supplied `benchmark_non_dyson.py`, real RHF references, spherical bases, identical geometry/frozen spaces, RHF tolerance 1e-12, and pole tolerance 1e-9. Processes ran serially without profiling. Eight threads were requested through the OMP/BLAS environment settings. NumPy uses Apple Accelerate; this local PySCF build lacks OpenMP. These are single-machine measurements, not universal speedups.

Process time below is wall time including interpreter startup and RHF. Setup time is the EPT constructor only, including integral transformation and static construction. A timeout supplies only a lower bound.

| Case | Original process (s) | Updated process (s) | Process speedup |
|---|---:|---:|---:|
| HF-cc-pvdz-ip-nD-NRL3 | 1.875 | 0.695 | 2.7x |
| N2-cc-pvdz-ip-nD-NRL3 | 17.866 | 1.329 | 13.4x |
| HF-cc-pvtz-ip-nD-NRL3 | >180 (stopped) | 8.238 | >21.9x |
| HF-cc-pvtz-ea-nD-NRL3 | 13.773 | 2.585 | 5.3x |

Updated EPT-only timings:

| Molecule / basis / sector / method | Setup (s) | Pole solves (s) | Status |
|---|---:|---:|---|
| HF-cc-pvdz-ip-nD-NRL3 | 0.267 | 0.005 | converged |
| N2-cc-pvdz-ip-nD-NRL3 | 0.991 | 0.014 | converged |
| HF-cc-pvtz-ip-nD-NRL3 | 7.854 | 0.014 | converged |
| HF-cc-pvtz-ip-NRL3 | 1.059 | 0.395 | converged |
| HF-cc-pvtz-ea-nD-NRL3 | 1.978 | 0.253 | converged |
| N2-cc-pvtz-ip-nD-NRL3 | 45.766 | 0.062 | converged |
| HF-cc-pvqz-ip-nD-NRL3 | — | — | timeout after 600.0 s |
| HF-aug-cc-pvtz-ip-nD-NRL3 | 167.002 | 0.036 | converged |

The smaller final nD-NRL3 eigenproblem does **not** guarantee an end-to-end speedup over NRL3. The reported NRL3 timing provides that comparison directly. Reuse the same `EPT` object for further targets at an unchanged reference so its static correction is not rebuilt.

## Numerical checks

The full test suite passes 132 cases. Checks include independent dense resolvents, indefinite sampling shifts, both spin sectors, normalized virtual-pair transformations, cached and streamed contractions, true-parent residual refinement, actual-residual failure rejection, derivatives, pole strengths, and AO normalization. Molecular diagnostics are JSON serializable.

Completed old/new comparisons:

| Case | Maximum energy change (eV) | Maximum PS change | Maximum static-matrix entry change (Ha) |
|---|---:|---:|---:|
| HF-cc-pvdz-ip-nD-NRL3 | 2.700e-13 | 3.697e-14 | 2.486e-12 |
| N2-cc-pvdz-ip-nD-NRL3 | 1.027e-12 | 9.492e-14 | 1.721e-11 |
| HF-cc-pvtz-ea-nD-NRL3 | 1.692e-13 | 1.110e-16 | 1.197e-12 |

A small linear-system residual does not certify good conditioning near an opposite-sector pole. Such cases can remain expensive and can raise `ConvergenceError`; increasing `static_max_cycle` is appropriate only when inspecting an iteration-limit failure. Broad timing guarantees across molecules, hardware, or basis sets are not claimed.

## Reproduce

Run with the same Python environment and resources against the old and updated source checkouts:

```sh
python examples/benchmark_non_dyson.py --basis cc-pvtz --threads 8 --output hf_tz.json
python examples/benchmark_non_dyson.py --basis aug-cc-pvtz --threads 8 --verbose 4 --output hf_aug_tz.json
python examples/benchmark_non_dyson.py --molecule N2 --basis cc-pvtz --threads 8 --output n2_tz.json
python -m pytest -q
```

The benchmark JSON includes inputs, environment versions, timing, pole energies/PS/residuals, and the full static matrix. The compact [recorded report](non_dyson_performance.json) omits the large static matrices but retains their old/new comparison errors and all completed pole data. The original source used for comparison is commit `39a7db4`; subsequent remote changes before this optimization concerned documentation/reference reports, not this solver.
