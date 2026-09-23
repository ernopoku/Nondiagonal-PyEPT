# Separate IP/EA representations

Two additional methods are available through `SectorEPT`. They solve separate
IP or EA eigenproblems and use transition matrices to obtain Dyson amplitudes.

| Name | Definition | Status |
|---|---|---|
| `nD-ADC(3)` | PySCF restricted non-Dyson IP/EA ADC(3), with full available transition moments | Established method; adapter checked against direct PySCF |
| `NRL3-ISR(3)` | Third-order canonical block reduction of this package's NRL3 auxiliary matrix | **Experimental construction**, derived below; dense small-system prototype |
| `nD-NRL3` | Symmetric static opposite-sector self-energy | Previous extension; remains under `EPT` |

**The experimental method is not a separately derived ground-state-metric
non-Dyson NRL3 theory.** It is a mathematically specified, NRL3-derived sector
representation. Neither its higher-order choices nor its equivalence to a
hypothetical direct NRL3 IP/EA metric derivation have external validation.
`nD-ADC(3)` is a different approximation from NRL3 and is never aliased to it.

## Installation and calculation

```bash
python -m pip install -e '.[test,sector]'
```

The ADC adapter requires **PySCF 2.14 or newer**. Existing EPT methods retain
their previous minimum dependency. Both new interfaces currently accept real,
canonical, molecular closed-shell RHF references and conventional integrals.

```python
from pyscf import gto, scf
from nondiagonal_ept import SectorEPT, run_sector_methods

mol = gto.M(
    atom="""
    F   0.0000   0.0000   0.0000
    H   0.0000   0.0000   0.9168
    """,
    unit="Angstrom", basis="6-31g", verbose=0,
)
mf = scf.RHF(mol).run(conv_tol=1e-12)

results = run_sector_methods(
    mf, ["nD-ADC(3)", "NRL3-ISR(3)"],
    frozen=1, sector="ip", nroots=3, tol=1e-9,
)
for method, poles in results.items():
    for root, pole in enumerate(poles):
        print(f"{method} root {root}: IP = {pole.binding_energy_ev:.6f} eV, "
              f"PS = {pole.strength:.6f}, residual = {pole.residual:.2e}")

attachment = SectorEPT(mf, "nD-ADC(3)", frozen=1, sector="ea")
for root, pole in enumerate(attachment.kernel(nroots=3)):
    print(f"root {root}: EA = {pole.binding_energy_ev:.6f} eV, PS = {pole.strength:.6f}")
```

Use `basis="cc-pvtz"` with `nD-ADC(3)` alone when appropriate for your system
and resources. The small basis above makes both methods suitable for a quick
example; it is not a basis-set recommendation for accurate affinities.

`nroots` counts states, **not original MO indices**. IP roots are ordered by
increasing ionization cost; EA roots by increasing attachment energy (largest
binding EA first). States may include satellites with small PS. `pole.target`
is `None`; `enumerate` supplies a zero-based root label. This API is separate
from `EPT.kernel(targets=[...])` and the existing JSON CLI. To compare a
specific orbital between methods, inspect Dyson character, not just root rank.
Degenerate orbitals may rotate within their degenerate subspace.

`frozen` uses the existing count or original spatial-MO-index-list convention.
`calc.original_mos` maps `pole.dyson_mo` back to original MO indices. The RHF
object is preserved. `spin=0/1` selects equivalent closed-shell spin channels.

## Energies, PS and visualization

Both methods return the existing `Pole` result type:

- `energy` is signed omega in hartree: minus IP for removal and minus binding
  EA for attachment. An unbound attachment has negative `binding_energy_ev`.
- `strength` is **one-spin PS**. PySCF's restricted ADC spectroscopic factor
  is spin-summed, so the adapter reports `P/2`. Its transition amplitudes
  already have the required normalization and are not divided by sqrt(2).
- `dyson_mo` and `dyson_ao` have squared MO/overlap-metric AO norm equal to PS.
  They include the transition matrix, rather than merely extracting the
  simple part of the sector eigenvector.
- `vector` is the sector eigenvector. ADC uses a spin-adapted representation;
  its Euclidean simple-component norm is not PS.
- `residual` is an independently evaluated eigenpair residual in hartree.
  ADC uses a relative Euclidean residual after its spin-adapted normalization.
  Failure raises `ConvergenceError`. `iterations=0` for ADC means the backend
  does not expose a count; the dense prototype reports `iterations=1`.

Cube export works as before:

```python
from pyscf.tools import cubegen
pole = results["nD-ADC(3)"][0]
cubegen.orbital(mol, "ip_root0.cube", pole.dyson_ao,
                nx=100, ny=100, nz=100, margin=5.0)
```

Display positive and negative isosurfaces. Do not square amplitudes before
exporting an orbital. Unit normalization is optional for comparing shapes,
not for reporting PS. See the Wiki's Dyson-orbital visualization guide.

## NRL3-specific derivation

This construction operates on the already defined **NRL3 auxiliary matrix**,
not directly on the physical many-electron Hamiltonian. Introduce a formal
coupling parameter lambda multiplying fluctuation integrals at fixed HF
orbital energies. Its exact polynomial decomposition is

```math
H(\lambda)=H_0+\lambda H_1+\lambda^2H_2+\lambda^3H_3.
```

Here H0 contains HF simple and triple orbital energies; H1 contains bare
simple–triple vertices and first-order triple interactions; H2 contains
NRL3's half-weight second-order vertices; H3 is the linear-singles static
simple block. The independent scaling test verifies this decomposition.
In particular, replacing these by ADC(3) blocks would change the method.

Partition the auxiliary space into

```math
P=\{\text{occupied simple},2h1p\},\qquad
Q=\{\text{virtual simple},2p1h\}.
```

Let D extract the two diagonal sector blocks and O the cross blocks. Seek
an antisymmetric generator, with zero PP and QQ blocks,

```math
S(\lambda)=\lambda S_1+\lambda^2S_2+\lambda^3S_3,
\qquad K=e^{-S}He^S.
```

Expanding the Baker–Campbell–Hausdorff series, define

```math
R_1=H_1,
```

```math
R_2=H_2+[H_1,S_1]+\tfrac12[[H_0,S_1],S_1],
```

```math
\begin{aligned}
R_3={}&H_3+[H_2,S_1]+[H_1,S_2]\\
 &+\tfrac12\left([[H_0,S_1],S_2]+[[H_0,S_2],S_1]+[[H_1,S_1],S_1]\right)\\
 &+\tfrac16[[[H_0,S_1],S_1],S_1].
\end{aligned}
```

At each order, cancel O(Kn) using the diagonal energies of H0:

```math
(S_n)_{pq}=-\frac{(R_n)_{pq}}{e_p^{(0)}-e_q^{(0)}}\quad(p\in P,q\in Q),
\qquad (S_n)_{qp}=-(S_n)_{pq}.
```

The effective sector matrices at lambda=1 are the PP and QQ blocks of

```math
K^{[3]}=H_0+D(R_1)+D(R_2)+D(R_3).
```

They are Hermitian and have no dynamical coupling to the other sector.
Cross-sector zero/near-zero denominators are rejected; neither clipping nor
energy shifting is used. Degeneracies *within* a sector need no such division.

### Consistent transition moments

Let J select all simple orbital rows of the original NRL3 auxiliary space.
Use an explicit **orthogonal completion** of the cubic generator,

```math
U=\exp(S_1+S_2+S_3),\qquad F_-=JUP,\quad F_+=JUQ.
```

For a normalized sector eigenvector y, the Dyson amplitudes are x=F y and
PS=x-transpose x. Through third order U agrees with

```math
I+S_1+(S_2+\tfrac12 S_1^2)
 +(S_3+\tfrac12(S_1S_2+S_2S_1)+\tfrac16 S_1^3).
```

Keeping the exponential adds specified fourth-and-higher-order terms to
transition moments. This is a deliberate completion, not conventional ADC's
transition-moment truncation. It guarantees nonnegative PS no larger than one
and the combined one-spin identity

```math
F_-F_-^T+F_+F_+^T=I.
```

The resulting approximate propagator is

```math
G^{\mathrm{ISR}}(z)=F_-(zI-K_-)^{-1}F_-^T+
                     F_+(zI-K_+)^{-1}F_+^T.
```

Away from singularities it agrees with J(zI-H)^(-1)J-transpose through third
order in lambda. It generally differs from fully resummed NRL3 starting at
fourth order. This is agreement with **NRL3**, not with an exact third-order
many-body expansion: NRL3 itself omits terms retained in strict ADC(3).
Combined completeness is not a proof of the correct separate electron-number
sum rules or of a physical correlated ground-state density.

### What remains unestablished

A direct non-Dyson NRL3 derivation from the article's Hermitized metric would
require separate N−1/N+1 metric and energy matrices, higher-order density
matrices, orthogonalization and transition operators, with a specified
consistent truncation. The canonical construction above does not supply that
independent many-electron derivation. The two could differ beyond the matched
orders. Size consistency of this new completed approximation and behavior
under strong correlation have not been established here.

The code is therefore explicitly named **NRL3-ISR(3)** and warns that it is
experimental. It is not aliased to `NRL3`, `nD-NRL3` or `nD-ADC(3)`.

## Cost and diagnostics

The ADC adapter uses PySCF's tensor-based sector matvecs and a residual-checked
Davidson solve, with the linear-dependence cutoff adjusted for tight requested
residuals. PySCF manages its integral storage. `max_memory_mb` is an estimate/
backend control, not an operating-system hard limit.

The NRL3 prototype constructs the **full auxiliary matrix first**, performs
dense commutators and a matrix exponential, and then solves the selected sector.
It currently offers **no speed or memory advantage** over matrix-free NRL3.
The default full-dimension guard is 1200, with an additional conservative
workspace estimate. Large cc-pVTZ calculations will usually exceed it. Raising
`max_dimension` requires enough `max_memory_mb` and can be expensive; prefer
small validation systems. `max_cycle/max_space` affect ADC only.

```python
calc = SectorEPT(mf, "NRL3-ISR(3)", frozen=1, sector="ip")
print(calc.diagnostics)
```

The diagnostics include generator norms, the minimum cross-sector denominator,
and norms of cross-sector coupling and higher-order terms left in the exactly
rotated parent matrix. These are truncation diagnostics, **not eigenpair
residuals or rigorous error bounds**. A tiny eigenpair residual only verifies
that the chosen approximate matrix was solved accurately.

## Validation

`tests/test_sector_isr.py` checks both ADC sectors against direct PySCF,
occupied/virtual freezing, energy signs, per-spin factors, full transition
amplitudes, AO norms, convergence failures, orbital phases and spin equality.
For NRL3-ISR(3), independent random-matrix exponentials verify fourth-order
remainders in the effective matrix, transition expansion and complete
propagator. Scaled molecular-integral tests verify fourth-order primary-pole
and Dyson-residue differences, plus the zero-interaction limit and combined
spectral completeness.

`python examples/validate_sector_methods.py` regenerates
[small-basis molecular comparisons](sector_comparison.json) for HF, water and
N2, both IP and EA. These compare implementations and approximations, not
accuracy against experiment. There is no independent published NRL3-ISR(3)
benchmark. See the current full-suite count in `test_results.txt`.

## References and provenance

- NRL3 blocks: Opoku, Pawłowski and Ortiz, *JCP* **159**, 124109 (2023),
  [doi:10.1063/5.0168779](https://doi.org/10.1063/5.0168779), especially Eq. 31
  and Table IV. That paper does **not** define the canonical prototype above.
- Non-Dyson ADC: Schirmer, Trofimov and Stelter, *JCP* **109**, 4734 (1998),
  [original article](https://files.isu.ru/ru/about/others/lab_kvant_himii/docs/nondys.pdf).
- [PySCF ADC documentation](https://pyscf.org/user/adc.html) and its cited
  implementation papers. The adapter uses PySCF; it does not independently
  reimplement or independently validate the established ADC tensor equations.
- General canonical block-diagonalization technique: Bravyi, DiVincenzo and
  Loss, *Annals of Physics* **326**, 2793–2826 (2011),
  [doi:10.1016/j.aop.2011.06.004](https://doi.org/10.1016/j.aop.2011.06.004).
  The NRL3 specialization and completion are the construction documented here,
  not a method claimed in that reference.
