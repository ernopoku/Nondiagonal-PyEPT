# BD-T1 fixed-core validation and correction

The supplied N2 reference uses a 1.1136 Angstrom bond, spherical cc-pVDZ,
two frozen occupied spatial orbitals, and unchanged RHF core orbitals throughout
Brueckner optimization. Python targets [6, 4, 3] correspond to the reference's
one-based active indices [5, 3, 2].

## Cause and correction

The previous BD-T1 method definition used full-weight linear doubles corrections
to the simple–triple coupling vertices. The supplied reference uses half weight:
B + C/2, rather than B + C. The four vertex weights are corrected. Triple-block
corrections and the Brueckner Fock block are unchanged. Other methods are unchanged.

The fixed-core orbital convention was already implemented by the PySCF BCCD
path. A new explicit invariant and regression test verify that frozen orbital
coefficients remain exactly equal to their input RHF values. No extra reference
keyword is needed: use EPT(mf, 'BD-T1', frozen=2) for this N2 example.

| Original MO (zero-based) | Supplied pole (Ha) | Before fix (Ha) | Corrected pole (Ha) |
|---|---:|---:|---:|
| 6 | -0.6139503411 | -0.6222124889 | -0.6139491792 |
| 4 | -0.5491151633 | -0.5711612561 | -0.5491134587 |
| 3 | -0.6742780381 | -0.7033247782 | -0.6742771185 |

Maximum pole difference after correction is 1.71e-6 Ha (about 0.000047 eV).
The converged pole strengths are 0.9293890, 0.9058776, and 0.8544687; the
supplied values are 0.9292955, 0.9059497, and 0.8545204. All Python eigenpair
residuals are below 1e-9 Ha. The BD total energy differs by approximately
2.22e-7 Ha; printed reference orbital energies also show small differences.
The input reference used less stringent stopping criteria. Exact equality of
these results is therefore not asserted, and the remaining differences have
not been individually apportioned between reference-orbital and solver errors.

Regression tolerances are 3e-6 Ha in pole energy, 2e-4 in pole strength, and
5e-7 Ha in BD total energy, alongside exact fixed-core and original-reference
preservation checks. Independent determinant-space tests verify the half-weight
linear vertex expansion. This benchmark covers N2 detachment, not a broad
BD-T1 molecular or attachment benchmark set.

## Reproduce

```bash
python examples/n2_bdt1.py > n2_bdt1_results.json
pytest -q
```

The recorded numerical report is [bdt1_n2_reference.json](bdt1_n2_reference.json).
Update an existing installation with `git pull` and `python -m pip install -e .`.
