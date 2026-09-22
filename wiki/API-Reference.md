# API reference

```python
EPT(mf, method="NRL3", frozen=0, sector="ip", spin=0,
    max_memory_mb=2000, brueckner_tol=1e-8)
calculation.kernel(targets=None, tol=1e-9, max_cycle=150, max_space=40)
calculation.dense_spectrum(max_dimension=2500)
calculation.self_energy(energy, derivative=False)
```


- `kernel(targets=[...])` follows the primary pole with largest overlap on each target simple operator. It returns only residual-converged roots or raises `ConvergenceError`. Strongly mixed primary/satellite states can require inspection of the full spectrum. Duplicate convergence from two targets raises an error.
- `dense_spectrum()` returns every pole for small systems, with a default dimension guard of 2500. It is also an independent eigensolver check on Davidson.
- `self_energy(E)` evaluates the non-diagonal Schur complement with checked MINRES solves. `self_energy(E, derivative=True)` returns its analytic energy derivative. Energies at triple-space singularities can fail to converge and raise an exception.
- A Dyson orbital is **not** unit-normalized by default. Its squared norm is its pole strength. `normalized_dyson_mo` provides the unit-normalized shape.
- The program uses spatial four-index ERIs and spin-block contractions, avoiding a full spin-orbital V⁴ tensor. It does not use density fitting or distributed/out-of-core propagator tensors. `max_memory_mb` is a conservative preflight estimate, not an operating-system memory limit. Very large basis sets remain expensive.


`targets=None` selects all active occupied spatial orbitals for IP and all active virtual spatial orbitals for EA. Explicit targets use the original zero-based spatial-MO indices. Both triple manifolds remain present.
