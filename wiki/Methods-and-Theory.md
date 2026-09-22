# Methods and theory

Both 2h1p and 2p1h operator manifolds are retained in a matrix-free Hermitian propagator Hamiltonian.

| Method | Approximation |
|---|---|
| ND2 / ADC(2) | Non-diagonal second-order Dyson self-energy |
| 2ph-TDA | First-order interactions in both triple manifolds |
| NR2 | Non-diagonal renormalized second order |
| NRP3 | Non-diagonal renormalized partial third order |
| NRQ3 | Non-diagonal renormalized quasiparticle third order |
| NRL3 | Non-diagonal renormalized linear third order |
| 3+ / ADC(3)-strict | Strict third-order Dyson ADC with ring/ladder renormalization |
| BD-T1 | Brueckner doubles with linear doubles/triple-operator terms |

[Full equations, tensor conventions, and block definitions](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/THEORY.md)

For asymmetric NR2/NRP3/NRQ3 truncations, `sector="ea"` selects the particle–hole counterpart. ND2, 2ph-TDA, NRL3, 3+, and BD-T1 share their Hamiltonian between IP and EA calculations.

The reference must be real molecular closed-shell RHF. BD-T1 prepares semicanonical Brueckner orbitals through PySCF BCCD. UHF, ROHF, DFT, density fitting, complex orbitals, periodic systems, and gradients are outside scope. Use `ADC(3)-strict`; the ambiguous name `ADC(3)` is rejected. Renormalized-static ADC(3) variants and diagonal-only methods are not included.

## Scientific references

- E. Opoku, F. Pawłowski, J. V. Ortiz, J. Chem. Phys. 159, 124109 (2023): [DOI](https://doi.org/10.1063/5.0168779).
- E. Opoku, F. Pawłowski, J. V. Ortiz, J. Chem. Phys. 155, 204107 (2021): [DOI](https://doi.org/10.1063/5.0070849).
- The [2021 article erratum](https://doi.org/10.1063/5.0167154) was identified, but its full text was unavailable during implementation. New-method definitions follow the supplied 2023 article and documented block definitions.

[[Home]] | [[Validation-and-Limitations]]
