# Citation

If you use Nondiagonal PyEPT in research, please cite the software and the
scientific articles relevant to the methods used. Citation is a scholarly
request, not an additional condition of the MIT license.

## Software

Ernest Opoku. *Nondiagonal PyEPT*, version 0.1.0. Python software.
https://github.com/ernopoku/Nondiagonal-PyEPT

Use GitHub's **Cite this repository** control for a formatted citation. Its
metadata come from [CITATION.cff](CITATION.cff). No software DOI has been assigned.
For reproducibility, also report the exact Git commit (`git rev-parse HEAD`),
PySCF version, geometry and units, basis, frozen orbitals, method, IP/EA sector,
targets, and convergence tolerance. Version 0.1.0 alone does not distinguish
updates made on the main branch.

```bibtex
@software{opoku_nondiagonal_pyept,
  author = {Opoku, Ernest},
  title = {Nondiagonal PyEPT: Reference Implementation of Non-diagonal Electron Propagator Methods},
  version = {0.1.0},
  url = {https://github.com/ernopoku/Nondiagonal-PyEPT}
}
```

## Scientific methods

- E. Opoku, F. Pawłowski, and J. V. Ortiz, *J. Chem. Phys.* **159**, 124109 (2023),
  [doi:10.1063/5.0168779](https://doi.org/10.1063/5.0168779).
- E. Opoku, F. Pawłowski, and J. V. Ortiz, *J. Chem. Phys.* **155**, 204107 (2021),
  [doi:10.1063/5.0070849](https://doi.org/10.1063/5.0070849).
  An erratum is available at [doi:10.1063/5.0167154](https://doi.org/10.1063/5.0167154).

See [theory](docs/THEORY.md) for implemented definitions and
[validation](docs/VALIDATION.md) for the scope of checks and remaining limitations.
Also acknowledge PySCF using its [citation guidance](https://pyscf.org/citation.html).


## Separate-sector methods

For `nD-ADC(3)`, also cite Schirmer, Trofimov and Stelter, *J. Chem. Phys.*
**109**, 4734 (1998), [article](https://files.isu.ru/ru/about/others/lab_kvant_himii/docs/nondys.pdf),
and the relevant implementation papers in the [PySCF ADC documentation](https://pyscf.org/user/adc.html).

For experimental `NRL3-ISR(3)`, cite the NRL3 paper and the exact software
commit, and identify it as the **canonical auxiliary-matrix construction**
defined in [SECTOR_ISR.md](docs/SECTOR_ISR.md). It is not a published
non-Dyson NRL3 theory. The general canonical transformation is described by
Bravyi, DiVincenzo and Loss, *Ann. Phys.* **326**, 2793 (2011),
[doi:10.1016/j.aop.2011.06.004](https://doi.org/10.1016/j.aop.2011.06.004).


## Conventional Dyson ADC(2) and ADC(3)

For these methods also cite J. Schirmer, L. S. Cederbaum and O. Walter,
*Phys. Rev. A* **28**, 1237 (1983),
[doi:10.1103/PhysRevA.28.1237](https://doi.org/10.1103/PhysRevA.28.1237).
For the `ADC(3)` DEM static correction, cite J. Schirmer and G. Angonoa,
*J. Chem. Phys.* **91**, 1754 (1989),
[doi:10.1063/1.457081](https://doi.org/10.1063/1.457081).
The equivalent resolvent/linear equations and particle-number limitation are
discussed by M. Deleuze, M. K. Scheller and L. S. Cederbaum,
*J. Chem. Phys.* **103**, 3578 (1995),
[doi:10.1063/1.470241](https://doi.org/10.1063/1.470241).
Report the explicit static scheme: `ADC(3)-DEM` is distinct from `3+` and
`nD-ADC(3)`. See [the definition](docs/DYSON_ADC.md).
