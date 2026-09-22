# Contributing to Nondiagonal PyEPT

Contributions to methods, numerical reliability, performance, examples, and
documentation are welcome. Please discuss substantial method or API changes in
an [issue](https://github.com/ernopoku/Nondiagonal-PyEPT/issues) first.

## Report a problem

Include a minimal runnable example, full traceback or unexpected output, package
versions, and the expected behavior. For numerical issues include Cartesian
geometry and units, basis, charge, frozen orbitals, method, IP/EA sector, target
indices, SCF tolerance, propagator tolerance, pole energies, and residual norms.
Remove private data from logs before posting. Explain the source and conventions
of comparison values; an energy match alone does not establish convergence.

## Development setup

Fork the repository, clone your fork, and create a branch for the change.
From the repository directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
pytest -q
```

Python 3.10 or newer is required. Use the virtual-environment activation command
appropriate to your shell and platform.

## Submit a pull request

1. Keep the change focused and describe the problem and resulting behavior.
2. Add independent numerical or algebraic checks for scientific changes. For a
   new method, document equations, tensor conventions, references, and its
   validation boundary in `docs/THEORY.md` and `docs/VALIDATION.md`.
3. Preserve zero-based original-MO target indexing, signed-pole conventions,
   residual-based convergence, and the existing single-method API unless an
   intentional interface change has been discussed.
4. Run relevant tests and report the commands, outcomes, and environment.
   Documentation-only changes need link and example checks, not numerical tests.
5. Update user-facing examples and documentation for interface changes. Wiki
   sources are mirrored in `wiki/`; maintainers synchronize changes to the live
   Wiki.

Do not silently replace reference data or relax tolerances to hide failures.
Distinguish algebraic tests, execution checks, and independent molecular
benchmarks when describing evidence.

## License and credit

Submit only material you have the right to contribute. Contributions to this
project are provided under its [MIT license](LICENSE). Keep applicable third-party
notices and cite scientific sources. Dependencies retain their own licenses.
Please communicate respectfully and give appropriate credit to collaborators.
