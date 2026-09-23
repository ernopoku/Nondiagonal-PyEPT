# API reference

```python
EPT(mf, method="NRL3", frozen=0, sector="ip", spin=0,
    max_memory_mb=2000, brueckner_tol=1e-8)
calculation.kernel(targets=None, tol=1e-9, max_cycle=150, max_space=40)
calculation.dense_spectrum(max_dimension=2500)
calculation.self_energy(energy, derivative=False)
```

## Solver behavior

`kernel` follows the primary pole with largest overlap on each target simple operator. It returns residual-converged roots or raises `ConvergenceError`. Duplicate convergence from two targets raises an error. Strong mixing may require inspecting the full spectrum.

`targets=None` selects active occupied spatial MOs for IP and active virtual spatial MOs for EA. Explicit targets are original zero-based spatial-MO indices. `spin=0` selects alpha and `spin=1` beta; their spectra agree for closed-shell references.

`dense_spectrum` includes satellites and is intended for small systems. Its default dimension guard is 2500.

`self_energy` evaluates the non-diagonal Schur complement through checked MINRES solves. With `derivative=True`, it returns the analytic energy derivative. Triple-space singularities can cause a convergence exception.

## Pole observables

- `energy`: signed propagator pole in hartree.
- `binding_energy_ev`: IP or EA, equal to minus the pole converted to eV.
- `strength`: squared norm of the Dyson orbital.
- `residual`: Hamiltonian eigenpair residual norm.
- `dyson_mo`: active spatial-MO coefficients.
- `dyson_ao`: AO coefficients, normalized in the AO overlap metric to the pole strength.
- `normalized_dyson_mo`: unit-normalized orbital shape.

## Memory

The implementation uses spatial four-index integrals and spin-block contractions. It avoids a full spin-orbital virtual fourth-order tensor, but remains an in-core research implementation. `max_memory_mb` is a conservative preflight estimate, not an operating-system memory cap.

[[Home]] | [[Examples]]
