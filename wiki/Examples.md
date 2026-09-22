# Examples

```python
from pyscf import gto, scf
from nondiagonal_ept import EPT

mol = gto.M(atom="H 0 0 0; F 0 0 0.9178", basis="cc-pvtz", verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)

calculation = EPT(mf, "NRL3", frozen=1, sector="ip")
for pole in calculation.kernel(targets=[4, 2], tol=1e-9):
    print(pole.binding_energy_ev, pole.strength, pole.residual)
    # pole.dyson_mo: active spatial-MO coefficients, norm squared = strength
    # pole.dyson_ao: AO coefficients, C.T @ S_AO @ C = strength
```

Targets are **zero-based indices in the ORIGINAL spatial MO list**, before freezing. For example, the fifth original spatial MO is always target `4`, including when a core orbital is frozen. `frozen=1` freezes the lowest occupied spatial MO; an explicit list can exclude occupied and/or virtual MOs. `spin=0` selects the alpha creator sector; `spin=1` selects beta. They have equal spectra for a closed-shell reference.

For attachment:

```python
attachment = EPT(mf, "NRQ3", frozen=1, sector="ea")
root = attachment.kernel(targets=[5])[0]
print(root.energy, root.binding_energy_ev)
```

`energy` is the signed propagator pole ω in hartree. For removal ω = E(N) − E(N−1), so IP = −ω. For addition ω = E(N+1) − E(N), so EA = −ω. Consequently, an unbound attachment has positive ω and negative EA. The program does not classify states using the sign of ω alone.

For BD-T1, use `EPT(mf, "BD-T1", frozen=1)`. This runs CCSD and orbital iterations, verifies the Brueckner singles norm, then builds the propagator. Target MO indices label the resulting semicanonical Brueckner orbitals; their character/order need not match the original canonical RHF orbitals. The supplied RHF object is preserved.


Command-line example from the repository root:

```bash
python -m nondiagonal_ept examples/hf.json -o hf_results.json
python examples/water.py
```
