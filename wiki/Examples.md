# Examples

## Ionization energies

```python
from pyscf import gto, scf
from nondiagonal_ept import EPT

mol = gto.M(atom="H 0 0 0; F 0 0 0.9178", basis="cc-pvtz", verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)
calculation = EPT(mf, "NRL3", frozen=1, sector="ip")
for pole in calculation.kernel(targets=[4, 2], tol=1e-9):
    print(pole.binding_energy_ev, pole.strength, pole.residual)
```

Targets use zero-based indices in the original spatial MO list, before freezing. `frozen=1` excludes the lowest occupied spatial MO; a list can exclude occupied and/or virtual orbitals.

## Electron attachment

```python
attachment = EPT(mf, "NRQ3", frozen=1, sector="ea")
root = attachment.kernel(targets=[5])[0]
print(root.energy, root.binding_energy_ev)
```

`energy` is the signed pole in hartree. For removal, omega = E(N) - E(N-1), and IP = -omega. For addition, omega = E(N+1) - E(N), and EA = -omega. An unbound attachment has positive omega and negative EA.

## Brueckner reference

```python
brueckner = EPT(mf, "BD-T1", frozen=1)
poles = brueckner.kernel()
```

BD-T1 runs orbital iterations and checks the singles norm. Targets then label semicanonical Brueckner orbitals, whose order and character may differ from RHF orbitals. The caller's RHF object is preserved.

## Command line

```bash
python examples/water.py
python -m nondiagonal_ept examples/hf.json -o hf_results.json
```

[[Home]] | [[API-Reference]]

## Cartesian geometry and multiple methods

Use one row per atom (element, x, y, z), with explicit ångström units. This example uses a 0.9168 Å bond, rather than the earlier 0.9178 Å bond.

```python
from pyscf import gto, scf
from nondiagonal_ept import run_methods

mol = gto.M(atom="""
F  0.0000  0.0000  0.0000
H  0.0000  0.0000  0.9168
""", unit="Angstrom", basis="cc-pvtz", verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)
results = run_methods(mf, ["NRL3", "NRQ3", "NRP3"],
                      frozen=1, sector="ip", targets=[4, 2], tol=1e-9)
for method, poles in results.items():
    print(method)
    for pole in poles:
        print(pole.target, pole.binding_energy_ev, pole.strength, pole.residual)
```

Update your checkout with `git pull` and reinstall with `python -m pip install -e .` before importing `run_methods`. It reuses the converged RHF reference, executes methods sequentially, and returns canonical method names mapped to lists of poles. Each method builds separate propagator intermediates. Duplicate methods are rejected; convergence failures raise an exception. All EPT constructor options and kernel solver controls are supported.

[Runnable example](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/examples/hf_multiple_methods.py)
