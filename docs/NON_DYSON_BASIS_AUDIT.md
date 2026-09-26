# Basis-set sensitivity of the static nD-NRL3 extension

**Historical report:** these measurements concern `full-space-v1`. The default is now the sector-projected version 2; see the [current method guide](NON_DYSON_NRL3.md) and [new validation](NON_DYSON_SECTOR_VALIDATION.md).

This investigation reproduces the reported discrepancy for HF at F (0,0,0), H (0,0,0.9168) Angstrom, RHF, one frozen occupied orbital, IP targets 4 and 2 (zero-based original MO indices). All calculations use the same reference for both methods. The source is commit `3696fc2`. RHF tolerance is 1e-12, pole tolerance 1e-10, and static tolerance 1e-11, tighter than the standard input.

## Reproduced results

| Basis | MO | NRL3 IP (eV) | nD-NRL3 IP (eV) | Difference (eV) |
|---|---:|---:|---:|---:|
| cc-pvdz | 4 | 15.248451 | 15.215615 | -0.032836 |
| cc-pvdz | 2 | 19.350583 | 19.340366 | -0.010217 |
| aug-cc-pvdz | 4 | 15.964862 | 15.926211 | -0.038651 |
| aug-cc-pvdz | 2 | 19.975676 | 20.017022 | +0.041347 |
| cc-pvtz | 4 | 15.822578 | 16.424596 | +0.602018 |
| cc-pvtz | 2 | 19.767618 | 20.217893 | +0.450276 |
| aug-cc-pvtz | 4 | 16.070782 | 16.341133 | +0.270351 |
| aug-cc-pvtz | 2 | 19.993542 | 22.776828 | +2.783286 |

The discrepancy persists with tighter tolerances. It is not explained by an unconverged eigenpair. The normalized Dyson-orbital overlaps and the full-parent static residuals are provided in [recorded data](non_dyson_basis_audit.json). Those checks test numerical solutions of the implemented equations; they do not validate the approximation against exact ionization energies.

## Mechanism in the implemented approximation

For IP, the code replaces the dynamic opposite-sector term

\[S^+(E)=B^+(EI-D^+)^{-1}(B^+)^T\]

by the real symmetric matrix

\[K_{pq}=\tfrac12\{S^+(\epsilon_p)_{pq}+S^+(\epsilon_q)_{pq}\}.\]

All active occupied and virtual simple orbitals remain in the matrix. Therefore an occupied-virtual entry K_ia includes evaluation at the **virtual** energy epsilon_a, even for a negative-energy IP root. This can sample a large opposite-sector resolvent response. In the spectral representation, denominators are epsilon_a minus an eigenvalue of D+; increasing the basis does not guarantee that these denominators remain well separated from zero. We have not computed the nearest D+ eigenvalue for every sampling point, so individual pole distances are not claimed here.

The numerical audit finds that the largest absolute occupied-virtual static entry grows substantially:

| Basis | max abs(K_ov), Ha | max abs(S+(eps_i)_ia), Ha |
|---|---:|---:|
| cc-pvdz | 0.063718 | 0.030817 |
| aug-cc-pvdz | 0.063601 | 0.032954 |
| cc-pvtz | 0.320594 | 0.029863 |
| aug-cc-pvtz | 1.021164 | 0.029011 |

The second column of values evaluates each occupied-virtual entry only at its occupied energy. It shows that the large entries are introduced by the virtual-energy sampling, rather than simply by an overall increase of the occupied-energy correction.

## Controlled diagnostic, not a replacement method

We changed only K_ia and its symmetric partner K_ai to S+(eps_i)_ia, leaving K_oo, K_vv, the retained 2hp block, couplings, and reference unchanged. This is an experimental matrix modification used to identify the source of the discrepancy. It is not the currently defined nD-NRL3 method, has not been derived as a complete non-Dyson theory, and is not enabled in production.

| Basis | MO | NRL3 IP (eV) | Diagnostic IP (eV) | Difference (eV) |
|---|---:|---:|---:|---:|
| cc-pvdz | 4 | 15.248451 | 15.230575 | -0.017875 |
| cc-pvdz | 2 | 19.350583 | 19.339116 | -0.011467 |
| aug-cc-pvdz | 4 | 15.964862 | 15.945226 | -0.019636 |
| aug-cc-pvdz | 2 | 19.975676 | 19.963092 | -0.012584 |
| cc-pvtz | 4 | 15.822578 | 15.805284 | -0.017294 |
| cc-pvtz | 2 | 19.767618 | 19.756332 | -0.011286 |
| aug-cc-pvtz | 4 | 16.070782 | 16.051423 | -0.019359 |
| aug-cc-pvtz | 2 | 19.993542 | 19.975516 | -0.018026 |

This controlled change strongly implicates the off-diagonal sampling prescription in these HF results. It does not establish a universally reliable corrected approximation; the virtual-virtual block and the EA sector still require a consistent treatment.

## Consequences

The double-zeta agreement is insufficient evidence of reliability across bases. The performance optimization preserves the chosen equations, including this limitation. More memory, more threads, or looser/tighter tolerances do not repair this approximation error. A future change of the sampling or sector-space definition must be explicitly named, theoretically justified, and validated across bases and molecules; it should not silently replace the current method.

For the reported triple-zeta cases, use the existing NRL3 calculation for the intended NRL3 predictions. Treat the current static nD-NRL3 results as exploratory. Manuscript claims of broad accuracy or basis-set robustness need to be revised to include this counterexample.

## Reproduction

Run `python examples/audit_non_dyson_basis.py` in the updated PyEPT environment. It writes per-basis JSON and NPZ matrices in `nd_nrl3_basis_audit_results` under the working directory. It uses internal Hamiltonian APIs for the controlled diagnostic and requires the performance update. The aug-cc-pVTZ case may take several minutes. No production method is modified.
