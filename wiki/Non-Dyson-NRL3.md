# Static non-Dyson NRL3

**Known basis-set limitation:** HF IP benchmarks show differences from NRL3 up to 0.60 eV with cc-pVTZ and 2.78 eV with aug-cc-pVTZ, despite converged residuals. The occupied–virtual static sampling is implicated. Double-zeta agreement does not establish general reliability. See the [basis-set investigation](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/NON_DYSON_BASIS_AUDIT.md) before using this exploratory approximation for quantitative predictions.

The method name is **`nD-NRL3`**. This project-defined extension freezes the opposite-sector NRL3 self-energy at HF orbital energies and symmetrizes its off-diagonal elements. It is distinct from a separately derived non-Dyson ADC intermediate-state formulation.

## Update your installation

```bash
git pull
python -m pip install -e .
```

## Run a calculation

```python
from pyscf import gto, scf
from nondiagonal_ept import EPT, run_methods

mol = gto.M(atom="""
F  0.0000  0.0000  0.0000
H  0.0000  0.0000  0.9168
""", unit="Angstrom", basis="cc-pvdz", verbose=0)
mf = scf.RHF(mol).run(conv_tol=1e-12)
calculation = EPT(mf, "nD-NRL3", frozen=1, sector="ip")
for pole in calculation.kernel(targets=[4, 2], tol=1e-9):
    print(f"MO {pole.target}: IP = {pole.binding_energy_ev:.6f} eV, "
          f"PS = {pole.strength:.6f}, residual = {pole.residual:.2e}")

results = run_methods(mf, ["NRL3", "nD-NRL3"], frozen=1,
                      sector="ip", targets=[4, 2], tol=1e-9)
```

For EA, set `sector="ea"` and choose virtual targets, e.g. `[5]` in this example. The JSON CLI also accepts `"method": "nD-NRL3"`. Existing original-MO indexing and frozen-core conventions apply.

## What is frozen?

For IP, 2ph contributions become static and only 2hp triples remain in the iterative Hamiltonian. For EA, reverse those roles. The full simple occupied and virtual orbital space is retained. For opposite-sector self-energy S, the static correction is

```text
K[p,q] = (S(eps_p)[p,q] + S(eps_q)[p,q]) / 2
```

This preserves Hermiticity and gives the diagonal frozen-at-own-HF-energy prescription. The opposite-sector contribution is retained, not discarded. Its linear systems are solved during setup with checked residuals. Setup still requires both sectors, so a smaller eigenproblem does not guarantee an overall speedup.

## PS and visualization

PS and `dyson_ao` are residues of this static approximate propagator. The frozen contribution has zero energy derivative, so PS generally changes relative to NRL3. No fitted intensity correction is used. Existing cube exports work unchanged. These residues are not the separately derived transition moments of non-Dyson ADC.

## Validation and limitations

The original validation included independent dense inverses, Hermiticity, dimension checks, analytic derivatives, solver-failure checks, spin equality, Dyson normalization, and the weak-coupling limit. Ten cc-pVDZ roots for HF, H2O, and N2 differ from NRL3 by 0.003715–0.142144 eV in absolute binding energy. The largest difference is the N2 target-6 IP. All nD-NRL3 eigenpair residuals are below 8.59e-10 Ha. The N2 IP secular dimension falls from 4016 to 761.

There is no external benchmark implementation for this extension. These are internal consistency checks and comparisons, not a claim of published-method validation. Near degeneracies, strong mixing, or opposite-sector resonances can cause larger changes. Small-basis attachment results are not basis-converged physical EA predictions.

[Complete equations and validation](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/NON_DYSON_NRL3.md) · [Numerical comparisons](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/non_dyson_comparison.json)

```bash
python examples/compare_non_dyson.py > comparison.json
pytest -q
```

[[Home]] | [[Methods and Theory|Methods-and-Theory]] | [[Visualizing Dyson Orbitals|Visualizing-Dyson-Orbitals]]

## Faster static setup

The optimized static solver uses preconditioned MINRES, evaluates only the needed
auxiliary sector, reuses packed nonzero spin and virtual-pair contractions, and stops residual refinement as soon
as the actual required residual is satisfied. The nD-NRL3
defaults are `static_tol=1e-10` and `static_max_cycle=5000` (per solve/refinement
pass). Set `mf.verbose = 4` before constructing `EPT` to see per-orbital progress.
`calculation.hamiltonian.static_diagnostics` records iterations, residuals,
operator products, cache size, and timing. The equations and PS normalization
are unchanged. Large or near-resonant static solves can still cost more than an
NRL3 calculation requesting only a few poles.

See the [performance report and benchmark script](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/NON_DYSON_PERFORMANCE.md).
