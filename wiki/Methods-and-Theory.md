# Methods and theory

Both 2h1p and 2p1h operator manifolds are retained in the matrix-free Hermitian Dyson propagator Hamiltonian. Separate-sector methods are described below.

| Method | Approximation |
|---|---|
| ND2 / ADC(2) | Non-diagonal second-order Dyson self-energy |
| ADC(3) | Conventional Dyson ADC(3) with Schirmer–Angonoa DEM static correction |
| 2ph-TDA | Two-particle-one-hole Tamm–Dancoff approximation |
| NR2 | Non-diagonal renormalized second order |
| NRP3 | Non-diagonal renormalized partial third order |
| NRQ3 | Non-diagonal renormalized quasiparticle third order |
| NRL3 | Non-diagonal renormalized linear third order |
| 3+ | Third-order plus |
| BD-T1 | Brueckner doubles with linear doubles/triple-operator terms |

[Full equations, tensor conventions, and block definitions](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/THEORY.md)

For asymmetric NR2/NRP3/NRQ3 truncations, `sector="ea"` selects the particle–hole counterpart. ND2, ADC(3), 2ph-TDA, NRL3, 3+, and BD-T1 share their Hamiltonian between IP and EA calculations.

The reference must be real molecular closed-shell RHF. BD-T1 prepares semicanonical Brueckner orbitals through PySCF BCCD. UHF, ROHF, DFT, density fitting, complex orbitals, periodic systems, and gradients are outside scope. `ADC(3)` now selects the DEM formulation; `ADC(3)-strict` remains an alias for `3+`. Full ADC(4), other static variants and diagonal-only methods are not included. See [[Dyson ADC|Dyson-ADC]] for examples and validation.

The static `nD-NRL3` extension is documented in [[Non-Dyson NRL3|Non-Dyson-NRL3]]. Established `nD-ADC(3)` and experimental `NRL3-ISR(3)` use `SectorEPT`; see [[Separate Sector Methods|Separate-Sector-Methods]].

## Scientific references

- E. Opoku, F. Pawłowski, J. V. Ortiz, J. Chem. Phys. 159, 124109 (2023): [DOI](https://doi.org/10.1063/5.0168779).
- E. Opoku, F. Pawłowski, J. V. Ortiz, J. Chem. Phys. 155, 204107 (2021): [DOI](https://doi.org/10.1063/5.0070849).
- J. Schirmer and G. Angonoa, J. Chem. Phys. 91, 1754 (1989): [DEM static self-energy](https://doi.org/10.1063/1.457081).
- The [2021 article erratum](https://doi.org/10.1063/5.0167154) was identified, but its full text was unavailable during implementation. New-method definitions follow the supplied 2023 article and documented block definitions.

[[Home]] | [[Validation-and-Limitations]]
