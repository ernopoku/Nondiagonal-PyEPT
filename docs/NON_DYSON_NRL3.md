# Static non-Dyson extension of NRL3

Invoke this implementation as **`nD-NRL3`**. It is the static opposite-sector
extension selected for this project. It is not a separately derived IP/EA ADC
intermediate-state representation, and it is not a new published benchmarked
method. The existing Dyson `NRL3` implementation is unchanged.

## Definition

Write the NRL3 self-energy as

\[
\Sigma(E)=\Sigma_\infty+B^- (EI-D^-)^{-1}(B^-)^T
                         +B^+ (EI-D^+)^{-1}(B^+)^T.
\]

The coupling matrices retain NRL3's half-weight linear corrections and the
triple blocks retain their first-order interactions. Let
`S_opp(E) = B_opp (E I - D_opp)^(-1) B_opp.T` and let eps_p be the canonical
HF orbital energy of simple orbital p. Freeze the opposite-sector contribution as

\[
K_{pq}=\tfrac12\{[S_{opp}(\epsilon_p)]_{pq}
                        +[S_{opp}(\epsilon_q)]_{pq}\}.
\]

For IP, the opposite sector is **2ph**; for EA, it is **2hp**. K is real symmetric.
Its diagonal is exactly S_opp(eps_p)[p,p], recovering the frozen-at-own-HF-energy
prescription on the diagonal. The static block is A_nd = A_NRL3 + K. The
remaining iterative Hamiltonian is

\[
H_{nd}=\begin{pmatrix} A_{nd}&B_{ret}\\B_{ret}^T&D_{ret}\end{pmatrix}.
\]

All active simple occupied and virtual orbitals of the selected spin remain.
Only the opposite triple manifold is removed from the secular problem. No
orbital-specific fitted energies or empirical parameters are introduced.
The static resolvent includes the fully renormalized opposite triple block,
not just its zeroth-order diagonal.

Each static solve uses diagonally preconditioned MINRES and checks the actual linear-system residual,
with residual refinement if needed. Failure raises ConvergenceError. A near
opposite-sector pole can make the chosen static approximation ill-conditioned;
small solver residuals alone do not establish physical reliability in that case.
The opposite blocks are needed during construction but not during the subsequent
Davidson iterations. Setup can therefore remain costly. No general speedup or
peak-memory reduction is claimed.

## Usage

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

For attachment use `sector="ea"` and suitable virtual targets (for this example,
`targets=[5]`). Targets retain the original zero-based spatial-MO indexing,
before freezing. The JSON CLI accepts `"method": "nD-NRL3"` with its existing
single-method input schema. The same RHF reference restrictions apply as for NRL3.

## Static setup controls and progress

The optimized implementation applies only the eliminated sector. It packs the
nonzero spin blocks of the 2ph contractions once and reuses BLAS-ready symmetric and antisymmetric
virtual-pair kernels. The latter cache is capped at the smaller of 1024 MB and
10% of `max_memory_mb`; larger tensors use streamed slabs instead. No dense
auxiliary Hamiltonian or spin-orbital four-virtual tensor is constructed.
The cache cap is additional-workspace management, not a hard operating-system
memory limit; the existing AO-to-MO memory estimate still applies.

The preconditioner is positive definite even for an indefinite shifted system:

```
gap = abs(eps_p - diagonal(D_opp))
floor = max(1e-8, 0.01 * max(gap))
M_inverse = 1 / maximum(gap, floor)
```

The floor regularizes only the **preconditioner**, not the self-energy denominator.
Periodic checks also stop MINRES once the actual residual meets the requested
threshold. This is especially important for tiny residual-refinement right-hand
sides, which the old fixed relative stopping rule could oversolve severely.
The final acceptance and residual-refinement checks use the original unpacked
parent operator, so even roundoff differences from pair compression are checked
and, if necessary, corrected. Every solution must satisfy:
`norm((eps_p I - D_opp) x - B_opp[p]) <= static_tol * max(1, norm(B_opp[p]))`.
Thus the optimization retains the same method, static matrix definition, and
residue normalization. A successful solver status with an inadequate residual
still raises an error. Conversely, reaching an iteration limit does not invalidate
a solution whose explicitly checked residual already meets the requested tolerance.

```python
mf.verbose = 4  # show per-orbital static setup progress
calculation = EPT(
    mf, "nD-NRL3", frozen=1, sector="ip",
    max_memory_mb=16000,
    static_tol=1e-10,
    static_max_cycle=5000,
)
poles = calculation.kernel(targets=[4, 2], tol=1e-9)
print(calculation.hamiltonian.static_diagnostics)
```

`static_tol` controls the shifted linear-system residuals; `kernel(tol=...)`
controls the final pole eigenpair residuals. `static_max_cycle` is the maximum
MINRES iteration count **per pass**; up to two residual-refinement passes are
allowed. Its nD-NRL3 default is 5000, independent of auxiliary dimension, replacing
the old dimension-dependent limit. Other methods keep their previous 500-cycle
default. Both controls also work through `run_methods` and JSON input.
Diagnostics include per-orbital iteration counts and residuals, total opposite-
block products, static setup/solve times, and ladder-cache bytes. With `verbose=0`,
setup remains quiet. Difficult near-pole shifts can still be expensive or fail;
a smaller final secular matrix does not guarantee a faster end-to-end calculation
than NRL3 when only a few poles are requested.

See [timings, validation, and reproduction commands](NON_DYSON_PERFORMANCE.md).

## PS, Dyson orbitals, and derivatives

PS and Dyson coefficients are the residues of this **static approximate
propagator**. They use the simple components of its normalized eigenvectors,
just as in the existing API. Thus `dyson_ao.T @ S @ dyson_ao == strength`, and
the existing cube-file visualization instructions apply.

The frozen K has zero energy derivative. `self_energy(E, derivative=True)`
therefore includes only the retained triple sector. PS values will generally
differ from full Dyson NRL3; no extra multiplicative correction is imposed.
These residues should not be identified with the separately derived transition
moments of non-Dyson ADC. The full simple space is retained, so the total spectral
weight over every reduced-Hamiltonian root is its simple-space dimension. This
is not a proof of a separate N-electron removal/addition sum rule.

## Checks and molecular comparisons

The original addition passed 55 tests; the optimized implementation and its additional checks are described in [the performance report](NON_DYSON_PERFORMANCE.md). New checks cover explicit
dense inversion of the opposite block for both sectors, off-diagonal
symmetrization, the diagonal freezing limit, reduced dimensions, Hermiticity,
analytic derivatives, the noninteracting limit, actual static-solve failures,
Davidson versus dense spectra, spin equality, PS and AO normalization, and
multi-method API operation. A nondegenerate weak-coupling check finds the
NRL3/nD-NRL3 energy difference decreases approximately 16-fold when the
fluctuation potential is halved, consistent with a fourth-order leading
change in that test. This does not establish exact third-order agreement with
full CI for NRL3 itself.

Ten IP/EA roots of HF, H2O, and N2 were compared at cc-pVDZ. The absolute
binding-energy differences from NRL3 range from **0.003715 to 0.142144 eV**.
The largest is the N2 target-6 IP. PS changes are also recorded; closeness of
energies does not imply identical intensities. All nD-NRL3 eigenpair residuals
are below 8.59e-10 Ha; static-solve residuals are below 9.45e-11 in these runs.
For N2, the matrix dimension falls from 4016 to 761 for IP and 3281 for EA.

These are reproducible internal consistency checks and comparisons, **not an
external validation of a published nD-NRL3 implementation**. The attachment
examples in this small basis should not be interpreted as converged predictions
of physical electron affinities. Diffuse basis convergence and broader benchmarks
remain to be studied. Near degeneracies, strong mixing, or opposite-sector
resonances can produce larger differences from NRL3.

Reproduce the comparisons with:

```bash
python examples/compare_non_dyson.py > comparison.json
pytest -q
```

[Recorded molecular data](non_dyson_comparison.json)

## Sources and scope of the extension

The diagonal non-Dyson prescription is described in Opoku, Pawlowski, and Ortiz,
J. Chem. Phys. 155, 204107 (2021), around Eq. (27),
[doi:10.1063/5.0070849](https://doi.org/10.1063/5.0070849). The NRL3 block definitions
come from J. Chem. Phys. 159, 124109 (2023),
[doi:10.1063/5.0168779](https://doi.org/10.1063/5.0168779). The symmetrized
non-diagonal static extension above is the explicit project definition; it is
not attributed to those articles as a published equation.
