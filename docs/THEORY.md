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
| BD-T1 | Brueckner Fock matrix | 1/2 | 1/2 | 1/2 | 1/2 | Linear in BD doubles | Linear in BD doubles |

For NR2/NRP3/NRQ3 attachment calculations, interchange holes and particles in this table. The supplied `rp3ea_tz.com` and `rq3ea_tz.com` use exactly that switch.

## BD-T1 triple corrections

The simple–triple vertices are `B + C/2`, with C as defined above. The triple-block corrections below are unchanged. Frozen orbitals remain fixed at the input RHF coefficients during optimization and semicanonicalization.

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

X_abcd is factorized through occupied pairs, not stored. Tests compare the formulas against the Hermitized linear-in-doubles contribution to the super-operator metric in small determinant spaces. A supplied N2/cc-pVDZ molecular reference is checked in the regression suite; see BD_T1_VALIDATION.md.

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

## Static opposite-sector extension

The project-defined `nD-NRL3` approximation now defaults to sector-projected version 2: occupied simple orbitals for IP, virtual simple orbitals for EA, and the corresponding retained triple sector. It freezes and symmetrizes the opposite-sector NRL3 resolvent only on that simple space. `static_space="full"` reproduces the legacy version. The full equations, energy-order argument, and residue convention are in [NON_DYSON_NRL3.md](NON_DYSON_NRL3.md).


## Conventional Dyson ADC additions

`ADC(2)` is identical to `ND2`. Standard `ADC(3)` retains the full-weight
second-order vertices and first-order triple blocks of `3+`, and replaces
its strict third-order static term by the Schirmer–Angonoa DEM static
self-energy. [DYSON_ADC.md](DYSON_ADC.md) defines the correlation-density
resolvents and linear response equation, including all spin factors and
validation. The `3+` implementation itself is unchanged.
