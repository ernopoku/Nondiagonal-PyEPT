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
  title = {Nondiagonal PyEPT},
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
