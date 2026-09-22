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
