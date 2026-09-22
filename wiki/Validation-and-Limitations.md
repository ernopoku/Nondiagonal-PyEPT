# Validation and limitations

All 38 automated tests passed for the publication copy. Tests cover independent determinant-operator algebra, PySCF MP2 agreement, the strict third-order limit against full CI, Davidson versus dense diagonalization, Hermiticity, spin symmetry, MO-phase invariance, spectral sum rules, Dyson normalization, and self-energy derivatives.

## Molecular comparisons

58 archived reference poles were compared for HF, H2O, N2, F2, and CO using spherical cc-pVTZ and the documented frozen cores. All computed residuals were below 1.895e-9 hartree. 53 energies agreed within 1e-5 hartree.

The remaining five differences were traced to premature energy-change-only stopping in the reference iterations. Reproducing the early iterations recovered the archived energies within about 1.3e-7 hartree, while eigenvector residuals remained large. Continuing to residual convergence changed the signed poles by up to approximately 0.068 eV. This is a convergence correction, not a measured improvement in physical accuracy.

[Detailed report, values, and reproduction commands](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/VALIDATION.md)

[Machine-readable comparisons](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/reference_comparison.json)

```bash
pytest -q
python examples/validate_reference.py --molecules HF H2O N2 F2 CO --output comparison.json
python examples/diagnose_convergence.py --molecule H2O --method NRL3 --target 3
```

## Limits

BD-T1 now has an N2/cc-pVDZ fixed-core molecular regression for three IP poles, agreeing within 1.8 microhartree. See [BD-T1 validation](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/BD_T1_VALIDATION.md). The articles' full statistical datasets have not been reproduced. The 2021 erratum full text was unavailable during implementation.

Only real molecular closed-shell RHF references are supported; BD-T1 prepares its own Brueckner reference. UHF/ROHF, DFT, density-fitted references, periodic systems, complex orbitals, gradients, non-Dyson ADC, diagonal-only methods, and renormalized-static ADC(3) variants are outside scope. Large-molecule and distributed performance has not been characterized.

[[Home]] | [[Methods-and-Theory]]
