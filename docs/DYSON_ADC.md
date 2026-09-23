# Conventional Dyson ADC(2) and ADC(3)

`EPT(mf, "ADC(2)")` computes the second-order Dyson self-energy, identical to
`ND2`. `EPT(mf, "ADC(3)")` computes standard Dyson ADC(3) with the
Schirmer–Angonoa **Dyson-expansion method (DEM)** for the static self-energy.
It is not an alias for `3+` or for non-Dyson ADC(3).

| Requested name | Formulation | Interface |
|---|---|---|
| `ADC(2)` / `ND2` | Second-order Dyson self-energy; zero static correlation term | `EPT` |
| `ADC(3)` / `ADC(3)-DEM` | Third-order dynamic self-energy with DEM static correction | `EPT` |
| `ADC(3)-strict` / `3+` | Same third-order dynamic blocks, strict third-order static correction | `EPT` |
| `nD-ADC(3)` | Separate IP/EA non-Dyson ADC(3), using PySCF | `SectorEPT` |

`Dyson-ADC(2)`, `ADC(2)-Dyson`, `Dyson-ADC(3)` and `ADC(3)-Dyson` are also
accepted. For backward compatibility the canonical name/result key for
`ADC(2)` remains `ND2`. Previously, the generic `ADC(3)` name was rejected;
it now explicitly selects the DEM formulation. Existing `3+` results are unchanged.

## Run one or several methods

```python
from pyscf import gto, scf
from nondiagonal_ept import EPT, run_methods

mol = gto.M(
    atom="""
    F   0.0000   0.0000   0.0000
    H   0.0000   0.0000   0.9168
    """,
    unit="Angstrom", basis="cc-pvtz", verbose=0,
)
mf = scf.RHF(mol).run(conv_tol=1e-12)

results = run_methods(
    mf, ["ADC(2)", "ADC(3)", "3+"],
    frozen=1, sector="ip", targets=[4, 2], tol=1e-9,
)
for method, poles in results.items():
    for pole in poles:
        print(f"{method} MO {pole.target}: IP = {pole.binding_energy_ev:.6f} eV, "
              f"PS = {pole.strength:.6f}, residual = {pole.residual:.2e}")

attachment = EPT(mf, "ADC(3)", frozen=1, sector="ea")
for pole in attachment.kernel(targets=[5]):
    print(f"EA = {pole.binding_energy_ev:.6f} eV, PS = {pole.strength:.6f}")
```

Targets remain zero-based **original spatial MO indices**, not energy-ranked
root indices. Both IP and EA are contained in the same Dyson Hamiltonian.
`sector` determines default targets and labels. Signed `pole.energy` is in
hartree; binding IP/EA is minus this pole. Unbound attachment has negative
binding EA. `pole.strength` remains the one-spin PS. MO/AO Dyson orbitals and
cube export work exactly as for the other `EPT` methods.

`python examples/dyson_adc.py` runs a smaller cc-pVDZ example. The JSON CLI
supports both methods and the static-solver controls:

```bash
python -m nondiagonal_ept examples/hf_adc3.json -o hf_adc3_results.json
```

## Static solver controls

```python
calc = EPT(mf, "ADC(3)", frozen=1,
           static_tol=1e-10, static_max_cycle=500)
print(calc.hamiltonian.static_diagnostics["static_residual"])
poles = calc.kernel(targets=[4, 2], tol=1e-9)
```

The static correction is computed once during construction. Its triple-space
resolvents and occupied–virtual response equation use matrix-free MINRES,
with an explicit residual check for every solve. `static_max_cycle` is the
iteration limit per static solve. Static controls are separate from
`kernel`'s `tol/max_cycle/max_space`. A failure raises `ConvergenceError`;
there is no silent fallback to `3+` or an unconverged static correction.

`static_diagnostics` includes residuals, the dynamic correlation-density
inhomogeneity and the completed DEM correlation density. The JSON CLI saves
scalar/list diagnostics, excluding the two density matrices.
`correlation_density_trace` is spin-summed. The DEM approximation need not
conserve particle number exactly; no trace adjustment or ad hoc rescaling
is applied. Its density is also not the density obtained by summing the final
Dyson removal poles. A residual is a numerical convergence test, not a bound
on physical errors.

## Working definition

The real RHF one-spin Dyson matrix has blocks

```math
\mathcal H=\begin{pmatrix}
\epsilon+\Sigma^\infty&B_h&B_p\\
B_h^T&D_h&0\\
B_p^T&0&D_p
\end{pmatrix}.
```

ADC(2) uses first-order vertices, zeroth-order triple energies and
Sigma-infinity = 0. ADC(3) uses full-weight second-order vertex corrections
and first-order triple interactions, exactly as `3+`. Only the static block
differs. The HF orbitals, energies and dynamic blocks are held fixed.

In DEM, retain `G0 + G0 (Sigma-infinity + M) G0` in the density contour
integral. Write the resulting correlation density as `Q + L(Sigma-infinity)`.
For occupied i,j and virtual a,b, define matrix-free resolvent solutions

```math
x_i^p=(\epsilon_i I-D_p)^{-1}B_{p,i}^T,\qquad
x_a^h=(\epsilon_a I-D_h)^{-1}B_{h,a}^T.
```

Then, with Q symmetric,

```math
Q_{ij}=-(x_i^p)^Tx_j^p,\qquad Q_{ab}=(x_a^h)^Tx_b^h,
```

```math
Q_{ia}=\frac{B_{p,a}x_i^p+B_{h,i}x_a^h}{\epsilon_i-\epsilon_a}.
```

The free-propagator static response has only mixed occupied–virtual blocks:

```math
[L(\Sigma)]_{ia}=[L(\Sigma)]_{ai}=\frac{\Sigma_{ia}}{\epsilon_i-\epsilon_a}.
```

For one-spin spatial density R and chemists' ERIs, let

```math
W(R)_{pq}=\sum_{rs}\big[2(pq|rs)-(pr|qs)\big]R_{rs}.
```

Solve the linear equation

```math
\Sigma^\infty=W\big(Q+L(\Sigma^\infty)\big).
```

The implementation solves only for `t_ia = Sigma_ia/(eps_i-eps_a)` using the
symmetric operator `gap*t - W(t+t-transpose)_ov`, then forms every static
block by a final W contraction. This is equivalent to the full spin-orbital
linear system, as independently checked in the tests. Frozen orbitals retain
their mean-field contribution; the DEM correlation/response space contains
only active orbitals.

Because the dynamic self-energy starts in second order and is complete
through third order, Q is correct through third order. Multiplication by W
and solution of the static response make the constant self-energy correct
through fourth order, with selected higher-order contributions. This does
**not** make the complete propagator ADC(4). In particular, this is not a
fully self-consistent iteration of the final Dyson density/Fock matrix and
is not the non-Dyson Sigma(4+) scheme.

## Validation and limitations

The added tests verify ADC(2)/ND2 identity; DEM Q against independent complex
contour integration; the static matrix against a full spin-orbital dense
linear solve; and its fourth-order accuracy against exact small-system FCI
correlation densities. Halving the fluctuation integrals produces fifth-order
static errors and fourth-order differences from `3+`. Molecular checks cover
IP/EA, frozen occupied and virtual MOs, alpha/beta equality, orbital phases,
Dyson norms, spectral completeness, the self-energy derivative residue
identity, and deliberate convergence failures.

[Reproducible molecular results](dyson_adc_comparison.json) include HF with
cc-pVDZ and cc-pVTZ, water/cc-pVDZ and N2/cc-pVDZ, using both IP and EA targets.
Generate them with `python examples/validate_dyson_adc.py`. These are internal
comparisons and regression evidence. No external program's standard
Dyson-ADC(3)-DEM molecular dataset has yet been compared; agreement with
non-Dyson ADC(3) is neither required nor used as a correctness criterion.

Real molecular canonical closed-shell RHF remains required. Zero HF gaps and
nonconvergent resolvent/response equations fail explicitly. There is no dense
full-Hamiltonian construction for this static correction; nevertheless,
large spatial integrals and triple tensors remain memory intensive.

## References

- J. Schirmer and G. Angonoa, *J. Chem. Phys.* **91**, 1754–1761 (1989),
  [doi:10.1063/1.457081](https://doi.org/10.1063/1.457081): DEM static self-energy.
- M. Deleuze, M. K. Scheller and L. S. Cederbaum, *J. Chem. Phys.* **103**,
  3578–3588 (1995), [doi:10.1063/1.470241](https://doi.org/10.1063/1.470241):
  linear static equations and ADC resolvents, Eqs. 9–10 and 38–41, and the
  particle-number limitation. The resolvent implementation above is written
  independently from these equations and checked by a separate contour oracle.
- J. Schirmer, L. S. Cederbaum and O. Walter, *Phys. Rev. A* **28**, 1237–1259
  (1983), [doi:10.1103/PhysRevA.28.1237](https://doi.org/10.1103/PhysRevA.28.1237).
- E. Opoku, F. Pawłowski and J. V. Ortiz, *J. Chem. Phys.* **159**, 124109
  (2023), [doi:10.1063/5.0168779](https://doi.org/10.1063/5.0168779):
  dynamic-block conventions and the distinction between ND2, strict `3+`
  and standard Dyson ADC(3).
