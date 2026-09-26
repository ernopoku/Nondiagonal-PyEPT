# Sector-projected static nD-NRL3 (version 2)

`EPT(mf, "nD-NRL3", ...)` now uses **sector-projected-v2** by default. This is an explicit revision of the original static extension, not a numerical retuning. IP keeps active occupied simple orbitals and 2hp configurations; EA keeps active virtual simple orbitals and 2ph configurations. The other triple sector contributes a static correction sampled only at the retained simple-orbital energies.

The original full-simple-space construction can be reproduced with `static_space="full"`. Its basis-set sensitivity is documented in the [legacy investigation](NON_DYSON_BASIS_AUDIT.md). Earlier numerical tables and manuscript dimensions for that formulation do not describe version 2. The existing Dyson NRL3 method is unchanged.

## Working equations

Write the parent NRL3 matrix, in a fixed spin sector, as

$$
H=\begin{pmatrix}A&B^-&B^+\\(B^-)^T&D^-&0\\(B^+)^T&0&D^+\end{pmatrix},
\qquad A=F+\Sigma_\infty.
$$

The vertices retain NRL3's half-weight linear corrections, and the triple blocks retain their first-order interactions. Let P select active occupied simple orbitals for IP, or active virtual simple orbitals for EA. R labels the retained triple sector (2hp for IP; 2ph for EA), and O labels the opposite triple sector. Define

$$
S_O(E)=B_O(EI-D_O)^{-1}B_O^T,
\qquad K_{pq}=\tfrac12\left([S_O(\epsilon_p)]_{pq}+[S_O(\epsilon_q)]_{pq}\right),\quad p,q\in P.
$$

The version-2 secular matrix is

$$
H_{\mathrm{nD},P}=\begin{pmatrix}
 A_{PP}+K_{PP}&B_{P,R}\\B_{P,R}^T&D_R
\end{pmatrix}.
$$

Implementation solves `(eps_p I - D_O) x_p = B_O[p,:].T` for each p in P and symmetrizes the resulting rows. No fitted shifts, denominator damping, adjustable energy cutoff, or molecule-specific branch is used. Occupation defines P; the sign of an orbital energy does not.

The matrix is the principal submatrix of the legacy static Hamiltonian on P plus the retained triple configurations. In particular, **K_PP is unchanged**. Version 2 removes the opposite-occupation simple orbitals; it does not clamp their energies or try to repair their large matrix elements. This also removes the virtual-energy sampling from IP static construction and occupied-energy sampling from EA static construction.

The entire active virtual basis still enters the integrals and triple configurations of an IP calculation. This is a projection of the simple-operator sector, not a frozen-virtual-orbital approximation. Frozen-core settings retain their original meaning.

At zeroth order this projection has a useful separation property. With highest active occupied energy epsilon_H, lowest active virtual energy epsilon_L, and a positive gap Delta = epsilon_L - epsilon_H,

$$
D^{+(0)}_{ia b}-\epsilon_p\ge 2\Delta\quad(p\ \mathrm{occupied}),
\qquad
\epsilon_p-D^{-(0)}_{ij a}\ge 2\Delta\quad(p\ \mathrm{virtual}).
$$

This follows from D+(0) = epsilon_a + epsilon_b - epsilon_i and D-(0) = epsilon_i + epsilon_j - epsilon_a. Thus the retained sampling energies lie on the appropriate side of the zeroth-order opposite-sector spectrum. The legacy full space did not have this property for its wrong-occupation rows. First-order interactions in D can still shift its spectrum; the inequality is not a guarantee for a strongly correlated, fully interacting resolvent.

## Perturbative justification and limits

For a canonical RHF reference with a finite occupied-virtual gap, scale the fluctuation potential by lambda. A selected-sector quasiparticle has a leading simple amplitude of order one, retained triple amplitudes of order lambda, and opposite-occupation simple amplitudes of order lambda squared: the simple-simple correlation coupling starts at second order, and the two-step path through triples also starts at second order. Their contribution to a nondegenerate target energy therefore starts at fourth order.

Also, replacing S_O(E) by S_O(eps_p) on the target diagonal changes the energy first at fourth order: the leading self-energy derivative is second order and E-eps_p is second order. Thus projection and freezing preserve the parent NRL3 target-energy expansion through third order under these assumptions. Independent weak-coupling tests verify fourth-order leading differences for both IP and EA.

This is an order argument about **energies**, not a claim that the full transition amplitudes, PS, or all matrix elements match NRL3 through third order. It is not an exact Schur-complement elimination of the removed simple space. That exact elimination would introduce additional energy-dependent blocks. Strong mixing, a closing gap, or a genuine opposite-sector resonance can invalidate the small-parameter argument.

This project-defined static construction is also distinct from the separately derived non-Dyson ADC intermediate-state representation. For that established approach use `SectorEPT(..., "nD-ADC(3)")`; see the [non-Dyson ADC literature](https://arxiv.org/abs/1910.07116). That literature is not a derivation or validation of this NRL3 extension.

## Use and compatibility

```python
from nondiagonal_ept import EPT

mf.verbose = 4
calculation = EPT(
    mf, "nD-NRL3", frozen=1, sector="ip",
    max_memory_mb=16000,
    static_space="sector",  # new default, explicitly shown for reproducibility
    static_tol=1e-10,
    static_max_cycle=5000,
)
for pole in calculation.kernel(targets=[4, 3, 2], tol=1e-9):
    print(pole.target, pole.binding_energy_ev, pole.strength)  # PS
print(calculation.hamiltonian.static_diagnostics)
```

Existing calls automatically use version 2 after updating the installation. Original zero-based MO indices are preserved. An IP target must be an active occupied orbital; an EA target must be an active virtual orbital. Frozen or wrong-sector targets raise an error.

`run_methods` and JSON input also accept `static_space`. Only nD-NRL3 supports a nondefault value. To reproduce the original approximation explicitly:

```python
legacy = EPT(mf, "nD-NRL3", frozen=1, sector="ip", static_space="full")
```

Diagnostics identify `formulation` (`sector-projected-v2` or `full-space-v1`), `static_space`, sampled original MO indices and energies, iteration counts, checked residuals, products, cache bytes, and timings. Progress output includes the formulation name. Save this metadata alongside results when comparing versions.

## Numerical controls and performance

`static_tol` controls the actual shifted-system residual; the default acceptance is

`norm((eps_p I - D_O) x_p - B_O[p,:]) <= 1e-10 * max(1, norm(B_O[p,:]))`.

Final acceptance is checked with the original uncompressed parent operator. Interactions in the unused, uncoupled triple sector on zero input are skipped during this check; the checked sector is unchanged. Positive-definite diagonal preconditioning, packed spin contractions, normalized virtual-pair kernels, and residual refinement remain enabled. `static_max_cycle` is a per-pass limit, with up to two refinement passes; it does not change the approximation. `kernel(tol=...)` independently controls the final eigenpair residual.

The ladder cache is capped at the smaller of 1024 MB and 10% of `max_memory_mb`; larger contractions use streamed slabs. The memory option is a working budget and integral-size check, not a hard process memory cap. No dense opposite-sector inverse is stored.

Only occupied shifts are needed for IP and only virtual shifts for EA. IP benefits particularly strongly because it no longer solves the high-energy virtual shifted systems. EA may still require many shifts. Integral transformations and vertex construction still involve the full active orbital space. Reuse one EPT object for all requested targets.

## PS, Dyson orbitals, and validation

`pole.strength` (printed as **PS**) is the squared norm of the retained simple component of a normalized eigenvector. `dyson_mo` retains the full active spatial-MO array shape, with zeros in the removed simple sector. `dyson_ao` is its AO expansion and satisfies `C.T @ S_AO @ C = PS`. These amplitudes belong to the projected approximate propagator; they are not the parent NRL3 amplitudes or separately corrected ADC transition moments.

`calculation.self_energy(E)` now returns a matrix in the retained simple space. Its row/column original-MO labels are in `static_diagnostics["original_mos"]`; do not index it directly by an original MO number.

Summing PS over the complete reduced spectrum gives the number of retained simple orbitals in the selected spin sector. This algebraic identity does not establish exact correlated IP/EA spectral sum rules.

See the [version-2 basis validation](NON_DYSON_SECTOR_VALIDATION.md) for molecular comparisons, residuals, PS, basis coverage, and failures if any. These comparisons use NRL3 as the parent approximation; agreement with it is not an experimental or exact-energy accuracy claim. No finite test set guarantees reliability for every molecule, basis, or strongly correlated reference.
