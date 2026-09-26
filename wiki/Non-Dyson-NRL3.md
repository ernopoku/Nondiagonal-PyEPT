# Non-Dyson NRL3: sector-projected version 2

Use `EPT(mf, "nD-NRL3", frozen=1, sector="ip")`. The default now keeps occupied simple orbitals and 2hp configurations for IP; EA keeps virtual simple orbitals and 2ph configurations. The opposite triple contribution is sampled only at the retained simple-orbital energies. This removes the wrong-sector sampling behind the original large triple-zeta discrepancies.

```python
calculation = EPT(mf, "nD-NRL3", frozen=1, sector="ip",
                  static_space="sector", max_memory_mb=16000)
poles = calculation.kernel(targets=[4, 3, 2], tol=1e-9)
print(calculation.hamiltonian.static_diagnostics["formulation"])
# sector-projected-v2
```

Original zero-based MO targets are unchanged. Existing scripts use version 2 after updating. To reproduce historical results, explicitly set `static_space="full"`; diagnostics then report `full-space-v1`. The revision changes the approximation, dimensions, energies, and PS; do not mix the two formulations in benchmark tables or manuscript claims.

Set `mf.verbose=4` for progress. The default static tolerance is 1e-10; `static_max_cycle=5000` is a per-pass iteration limit. All targets on one EPT object share the same static setup. Dyson MO arrays retain their active-space length, with zeros outside the retained simple sector. Printed pole strengths are PS.

See the [complete equations and theoretical limits](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/NON_DYSON_NRL3.md), [basis validation and reproduction instructions](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/NON_DYSON_SECTOR_VALIDATION.md), and [historical basis-set diagnosis](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/NON_DYSON_BASIS_AUDIT.md).

This is a static project-defined extension of NRL3, not the separate non-Dyson ADC ISR. Validation against the parent method does not guarantee experimental accuracy or reliability for every basis and molecule.
