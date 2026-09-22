# Methods and theory


| Name | Meaning |
|---|---|
| `ND2` / `ADC(2)` | Non-diagonal second-order Dyson self-energy |
| `2ph-TDA` | First-order interactions in both triple manifolds |
| `NR2` | Non-diagonal renormalized second order |
| `NRP3` | Non-diagonal renormalized partial third order |
| `NRQ3` | Non-diagonal renormalized quasiparticle third order |
| `NRL3` | Non-diagonal renormalized linear third order |
| `3+` / `ADC(3)-strict` | Strict third-order Dyson ADC with ring/ladder renormalization |
| `BD-T1` | Brueckner-doubles reference with terms linear in doubles and triple operators |

`sector="ea"` applies the particle–hole counterpart of the asymmetric NR2/NRP3/NRQ3 truncations. ND2, 2ph-TDA, NRL3, 3+, and BD-T1 have the same Hamiltonian for IP and EA. Merely changing the sign of an IP does not implement an EA-specific NRQ3 calculation.

**Reference scope:** real, molecular, closed-shell RHF; BD-T1 automatically constructs a semicanonical Brueckner reference with PySCF BCCD. UHF, ROHF, DFT, density-fitted SCF references, complex/spinor orbitals, periodic systems, gradients, and non-Dyson ADC are not implemented. The generic name `ADC(3)` is deliberately rejected: use `ADC(3)-strict`. Fourth-order/renormalized-static ADC(3) variants and the articles' diagonal-only methods are outside this implementation.


# Equations and conventions

The 2023 article's Table IV specifies the new non-diagonal methods. This document records the equations actually evaluated; it is not a transcription of the articles' diagonal-method catalog.

## Index and normalization conventions

Real orthonormal spin orbitals are ordered as consecutive alpha/beta pairs of each active spatial orbital, with all occupied pairs before virtual pairs. Indices i,j,k,l are occupied; a,b,c,d are virtual; p,q,r,s span all active spin orbitals.

Let `g[p,q,r,s] = <pq||rs>` in antisymmetrized Dirac notation:

\[
g_{pqrs}=(pr|qs)\delta_{\sigma_p\sigma_r}\delta_{\sigma_q\sigma_s}
 -(ps|qr)\delta_{\sigma_p\sigma_s}\delta_{\sigma_q\sigma_r}.
\]

The propagator basis consists of simple creators, 2hp creators `a_i† a_j† a_a` with i<j, and 2ph creators `a_a† a_b† a_i` with a<b. Each distinct triple has unit zeroth-order metric. No extra √2 belongs on these unique-pair matrix elements. Tensor actions temporarily unpack the unique pairs antisymmetrically; sums over all ordered pairs consequently carry 1/2.

Only operators with net alpha creation, or net beta creation when requested, enter a calculation. RHF alpha and beta spectra are degenerate. A spin-restricted reference does not justify dropping opposite-spin intermediate configurations.

## Reference amplitudes

The MP partition uses canonical energies ε and the fluctuation potential. With the excitation-operator convention of 2023 Eqs. 32–33,

\[
t_{ij}^{ab}=g_{ijab}/(\epsilon_i+\epsilon_j-\epsilon_a-\epsilon_b),
\]

\[
s_i^a=\frac{1}{2(\epsilon_i-\epsilon_a)}
\left[\sum_{jkb}t_{jk}^{ab}g_{bijk}
+\sum_{jbc}t_{ij}^{bc}g_{bcaj}\right].
\]

`t` is first order and `s` second order. Near-zero denominators are rejected rather than shifted silently. The spin MP2 energy `1/4 sum(g_ijab*t_ijab)` is checked against PySCF MP2.

## Simple–triple blocks

The first-order couplings are

\[
B^-_{p,ija}=g_{paij},\qquad B^+_{p,iab}=g_{piab}.
\]

Write the full-weight second-order vertex pieces as

\[
C^-_{p,ija}=\tfrac12\sum_{bc}g_{pabc}t_{ij}^{bc}
+\sum_{kb}(g_{pkib}t_{jk}^{ab}-g_{pkjb}t_{ik}^{ab}),
\]

\[
C^+_{p,iab}=\tfrac12\sum_{kl}g_{pikl}t_{kl}^{ab}
+\sum_{kc}(g_{pcak}t_{ik}^{bc}-g_{pcbk}t_{ik}^{ac}).
\]

The old symmetric metric uses `B+C`. The new Hermitized, intermediately normalized metric uses `B+C/2`. The first term in each C is the ladder contraction; the other two are the antisymmetrized ring contractions.

Independent tests evaluate the same quantities as reference/first-order-double matrix elements of `{[a_p,H], triple}` in a finite determinant space, without using these tensor formulas.

## Simple–simple block

\[
A_{pq}=\epsilon_p\delta_{pq}+\Sigma^\infty_{pq},\qquad
\Sigma^\infty_{pq}=\sum_{rs}g_{prqs}\,\Delta\rho_{rs}.
\]

For NRP3/NR2/ND2/2ph-TDA, Δρ=0. For NRQ3/NRL3, only `Δρ_ia=Δρ_ai=s_i^a/2` is retained. For 3+,

\[
\Delta\rho_{ij}=-\tfrac12\sum_{kab}t_{ik}^{ab}t_{jk}^{ab},\quad
\Delta\rho_{ab}=\tfrac12\sum_{ijc}t_{ij}^{ac}t_{ij}^{bc},\quad
\Delta\rho_{ia}=\Delta\rho_{ai}=s_i^a.
\]

The diagonal double-amplitude terms have zero total trace together. They are deliberately absent from the new linear metric. Static density contractions are independently checked against explicit reference-state operator matrix elements.

BD-T1 instead uses the determinant Fock matrix of the converged Brueckner orbitals. Its oo and vv blocks are semicanonical; its ov block need not vanish. MP2 singles or a quadratic MP2 density must not be substituted for that block.

## Triple–triple blocks

For antisymmetric tensors x_ija and y_iab,

\[
(D^-_0x)_{ija}=(\epsilon_i+\epsilon_j-\epsilon_a)x_{ija},
\quad
(D^+_0y)_{iab}=(\epsilon_a+\epsilon_b-\epsilon_i)y_{iab}.
\]

The first-order interactions act as

\[
(V^-x)_{ija}=-\tfrac12\sum_{kl}g_{ijkl}x_{kla}
+\sum_{kb}(g_{biak}x_{kjb}-g_{bjak}x_{kib}),
\]

\[
(V^+y)_{iab}=\tfrac12\sum_{cd}g_{abcd}y_{icd}
-\sum_{jc}(g_{jaic}y_{jcb}-g_{jbic}y_{jca}).
\]

These actions are checked against projected physical Hamiltonians in N−1 and N+1 determinant spaces, including the appropriate E_HF subtraction and IP sign. The four-virtual contraction is performed in spatial spin blocks: antisymmetry makes its exchange part equal to its direct part after the factor 1/2.

No 2hp–2ph interaction block is included in the named methods.

## Method switches for ionization

Each vertex cell gives the multiplier of C added to B. D entries give the highest interaction order retained.

| Method | Static density | h–2hp | p–2hp | h–2ph | p–2ph | D− | D+ |
|---|---|---:|---:|---:|---:|---:|---:|
| ND2 | None | 0 | 0 | 0 | 0 | 0 | 0 |
| 2ph-TDA | None | 0 | 0 | 0 | 0 | 1 | 1 |
| NR2 | None | 1/2 | 0 | 0 | 0 | 1 | 0 |
| NRP3 | None | 1/2 | 1/2 | 0 | 0 | 1 | 0 |
| NRQ3 | Linear singles | 1/2 | 1/2 | 0 | 0 | 1 | 0 |
| NRL3 | Linear singles | 1/2 | 1/2 | 1/2 | 1/2 | 1 | 1 |
| 3+ | Singles + quadratic doubles | 1 | 1 | 1 | 1 | 1 | 1 |
| BD-T1 | Brueckner Fock matrix | 1 | 1 | 1 | 1 | Linear in BD doubles | Linear in BD doubles |

For NR2/NRP3/NRQ3 attachment calculations, interchange holes and particles in this table. The supplied `rp3ea_tz.com` and `rq3ea_tz.com` use exactly that switch.

## BD-T1 triple corrections

In this section t denotes the **converged Brueckner doubles**, not MP2 t. Define

\[
r_{ij}=\tfrac12\sum_{kab}(g_{ikab}t_{jk}^{ab}+g_{jkab}t_{ik}^{ab}),
\]
\[
r_{ab}=\tfrac12\sum_{ijc}(g_{ijac}t_{ij}^{bc}+g_{ijbc}t_{ij}^{ac}),
\]
\[
X_{ijkl}=\tfrac12\sum_{ab}(g_{ijab}t_{kl}^{ab}+g_{klab}t_{ij}^{ab}),
\]
\[
X_{abcd}=\tfrac12\sum_{ij}(g_{ijab}t_{ij}^{cd}+g_{ijcd}t_{ij}^{ab}),
\]
\[
W_{iajb}=\sum_{kc}(g_{ikac}t_{jk}^{bc}+g_{jkbc}t_{ik}^{ac}).
\]

The added triple actions are

\[
(\Delta D^-x)_{ija}=-\tfrac14\sum_{kl}X_{ijkl}x_{kla}
-\tfrac12\sum_{kb}(W_{iakb}x_{kjb}-W_{jakb}x_{kib})
+\tfrac12\sum_l(x_{ila}r_{jl}-x_{jla}r_{il})
+\tfrac12\sum_bx_{ijb}r_{ab},
\]

\[
(\Delta D^+y)_{iab}=\tfrac14\sum_{cd}X_{abcd}y_{icd}
+\tfrac12\sum_{jc}(W_{iajc}y_{jcb}-W_{ibjc}y_{jca})
-\tfrac12\sum_c(y_{icb}r_{ac}-y_{ica}r_{bc})
-\tfrac12\sum_jy_{jab}r_{ij}.
\]

X_abcd is factorized through occupied pairs, not stored. Tests compare the formulas against the Hermitized linear-in-doubles contribution to the super-operator metric in small determinant spaces. Independent published BD-T1 molecular benchmark validation remains outstanding.

## Poles, self-energy, and intensity

The enlarged Hamiltonian is

\[
\mathcal H=\begin{pmatrix}
A&B^-&B^+\\(B^-)^T&D^-&0\\(B^+)^T&0&D^+
\end{pmatrix}.
\]

Diagonalizing it is equivalent to solving the non-diagonal Dyson equation with

\[
\Sigma(E)=A-\epsilon+B(EI-D)^{-1}B^T.
\]

The full inverse is retained, including higher-order terms. In particular, `(B+C/2) R (B+C/2)^T` is not subsequently truncated to remove its quadratic-vertex contribution; this follows the 2023 non-diagonal Hamiltonian prescription. It should not be confused with directly implementing a truncated or diagonal-only formula from the 2021 article.

A unit-normalized enlarged eigenvector has simple part c. The residue/pole strength is `Z=c.T@c`, and the unnormalized Dyson orbital is `sum_p c_p phi_p`. With u=c/√Z,

\[
\Sigma'(E)=-B(EI-D)^{-2}B^T,\qquad
Z=[1-u^T\Sigma'(\omega)u]^{-1}.
\]

The implementation verifies the AO overlap norm and tests the derivative against finite differences. For the full finite spectrum, `sum_n Z_n` equals the number of simple operators in the selected spin sector. This sum rule does not imply that each approximate pole is experimentally accurate.

The production convergence condition is `||H C − ω C||₂ < tol`. A small change in successive Ritz energies alone is insufficient, as demonstrated by the archived water and nitrogen calculations in the validation report.

## References

- E. Opoku, F. Pawłowski, J. V. Ortiz, *J. Chem. Phys.* **159**, 124109 (2023), [doi:10.1063/5.0168779](https://doi.org/10.1063/5.0168779), especially Eqs. 18–24, 28–33 and Table IV.
- E. Opoku, F. Pawłowski, J. V. Ortiz, *J. Chem. Phys.* **155**, 204107 (2021), [doi:10.1063/5.0070849](https://doi.org/10.1063/5.0070849).
- An erratum exists at [doi:10.1063/5.0167154](https://doi.org/10.1063/5.0167154). Its full text was not available during this work; the implemented new-method definitions use the supplied 2023 article and explicitly documented block definitions, rather than assuming the uncorrected 2021 formulas are definitive.
- [PySCF AO-to-MO documentation](https://pyscf.org/contributor/ao2mo_developer.html). The BD reference uses `pyscf.cc.bccd.bccd_kernel_` and the spin-amplitude conversion in `pyscf.cc.addons`.
