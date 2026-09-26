# 0.2.0: explicit nD-NRL3 formulation revision

`nD-NRL3` now defaults to the sector-projected version 2 described in [the method guide](NON_DYSON_NRL3.md). This intentionally changes numerical results, Hamiltonian dimensions, and Dyson amplitudes. It corrects the wrong-sector static sampling responsible for the legacy triple-zeta discrepancies; it is not an exact reformulation of the full-space approximation.

- IP retains occupied simple orbitals plus 2hp configurations; EA retains virtual simple orbitals plus 2ph configurations.
- Opposite-sector static solves sample only retained simple-orbital energies.
- Original MO target numbering, frozen-core definitions, resource settings, and tolerances are unchanged.
- Wrong-sector targets are rejected. Dyson MO arrays retain their full active-orbital length, with zeros outside the selected simple sector.
- `static_space="full"` explicitly reproduces version 1. The new default is `static_space="sector"`; Python, `run_methods`, and JSON input support it.
- Saved static diagnostics and verbose progress identify the formulation. Preserve this metadata for reproducibility.
- Other propagator methods are unchanged.

Update the same Python environment used to run calculations:

```sh
cd ~/Nondiagonal-PyEPT
git pull --ff-only
python -m pip install -e .
```

Existing nD-NRL3 scripts then use version 2. For manuscript work, update the working equations and recompute result and dimension tables; old full-space tables should not be relabeled as version 2. See [validation and its limits](NON_DYSON_SECTOR_VALIDATION.md).
