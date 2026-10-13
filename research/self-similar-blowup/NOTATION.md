# Shared conventions for the self-similar blowup project

All agents write LaTeX *fragments* (no preamble, no `\begin{document}`) into
`sections/`, code into `code/<prefix>_*.py`, figures (PDF or PNG) into
`figures/<prefix>_*.pdf`, BibTeX into `bib/<prefix>.bib`.
Do **not** edit files owned by other agents. Do **not** run git.

The master document loads: `amsmath, amssymb, amsthm, mathtools, bm,
graphicx, booktabs, hyperref, cleveref, enumitem`. Theorem environments
available: `theorem, proposition, lemma, corollary, definition, remark,
conjecture, example` (all numbered within section). Use `\label{<prefix>:...}`.
Each fragment starts with `\section{...}`; use `\subsection` below that.
Figures: `\includegraphics[width=...]{figures/<file>}` (path relative to the
project root `research/self-similar-blowup/`).

BibTeX keys: `<Lastname1><Lastname2>` or `<Lastname1>EtAl` followed by year,
ASCII only, e.g. `DonningerSchoerkhuber2012`, `CostinDonningerGlogic2017`,
`Donninger2011`, `GlogicSchoerkhuber2021`, `MerleZaag2007`. Add `a/b` for
same authors same year. Only cite papers you have verified exist (title,
authors, journal/arXiv, year). If unsure, write the claim without a citation
and flag it with `% TODO verify` in the source.

## Coordinates (blowup point (T, x0 = 0))

Equation  `u_tt - Δu = F(u)` on `R^{1+d}` (lower-order terms such as `γ u_t`,
`m^2 u` are treated as perturbations, see below).

Similarity coordinates on the backward light cone
`Γ_T = {(t,x): 0 ≤ t < T, |x| ≤ T - t}`:

    τ = -log(T - t) + log T,      ξ = x / (T - t),   so (τ,ξ) ∈ [0,∞) × \bar B^d.

`∂_τ|_ξ = (T-t)∂_t - x·∇_x`,  `∇_ξ = (T-t)∇_x`,  `Λ := ξ·∇_ξ`.

### Power-type scaling (exponent s)
If `F` is homogeneous so that `u ↦ μ^{s} u(μ t, μ x)` is a symmetry
(power: `s = 2/(p-1)`; MEMS `F = -λ w^{-q}`: `s = -σ`, `σ = 2/(q+1)`), set

    ψ1 = (T-t)^{s}   u,         ψ2 = (T-t)^{s+1} ∂_t u.

Then (check this yourself)

    ∂_τ ψ1 = -Λψ1 - s ψ1 + ψ2,
    ∂_τ ψ2 =  Δψ1 - Λψ2 - (s+1) ψ2 + (T-t)^{s+2} F(u),

and `(T-t)^{s+2}F(u) = F(ψ1)` for the homogeneous nonlinearities above.

### Exponential scaling
For `F = δ e^u`: `ψ1 = u + 2 log(T-t)`, `ψ2 = (T-t) ∂_t u`;
`∂_τψ1 = -Λψ1 + ψ2 - 2`, `∂_τψ2 = Δψ1 - Λψ2 - ψ2 + δ e^{ψ1}`.

### Lower-order terms
`γ u_t` contributes `-γ e^{-τ}` × (ψ2-type term) (rate e^{-τ});
`m^2 u` contributes `-m^2 e^{-2τ} ψ1` (rate e^{-2τ}), with `T=1`
normalisation `T - t = e^{-τ}`.

## Linearisation around the ODE (spatially homogeneous) profile

If `u_T(t) = κ (T-t)^{-s}` (or the log analogue) and `V0 := F'(u_T)(T-t)^2`
(a constant), eigenmodes `(ψ1,ψ2) = e^{λτ}(g, ...)` of the linearised
first-order operator correspond to `v = (T-t)^{-s-λ} g(ξ)` solving
`v_tt - Δv = V0 (T-t)^{-2} v`, i.e. with `a := s + λ`,

    (δ_ij - ξ_i ξ_j) ∂_i∂_j g - 2(a+1) Λ g - (a(a+1) - V0) g = 0   in B^d.

(For the exponential case use s = 0, V0 = 2.) Writing `g = r^ℓ Y_ℓ(ω) h(r^2)`
gives a hypergeometric equation, `h = 2F1(α,β; d/2+ℓ; r^2)` with
`α = (a+ℓ-a1)/2, β = (a+ℓ-a2)/2`, where `a1,a2` are the roots of
`a^2 + a - V0 = 0`. Smoothness at r = 1 forces α ∈ -N0 or β ∈ -N0 (modulo
the integer/log cases, which must be checked), giving eigenvalues

    λ = a1 - s - ℓ - 2n,   λ = a2 - s - ℓ - 2n,   n ∈ N0.

This is a *hint to be verified independently*, not a result to be copied.
Pay attention to: the function space (H^k × H^{k-1}(B^d), required k), the
growth bound of the free operator, the log/integer cases, the symmetry modes
(time translation → λ = 1; Lorentz boosts → λ = 0, ℓ = 1; other symmetries),
and whether eigenvalues are semisimple.

## Donninger et al. template (what each model section should contain)
1. Model, physical/industrial origin, nondimensionalisation.
2. Self-similar (ODE and, if any, non-trivial) blowup profiles; symmetry
   families (time translation, space translation, Lorentz boosts, scaling,
   shift).
3. First-order system in similarity coordinates; linearised operator
   `L = L0 + L'`; semigroup generation and growth bound for `L0` on
   `H^k × H^{k-1}(B^d)` (Lumer–Phillips with a suitable inner product, cf.
   Donninger–Schörkhuber, Ostermann, Glogić); compactness of `L'`.
4. Spectral analysis: unstable/neutral eigenvalues, their symmetry origin,
   mode stability (explicit hypergeometric/connection argument or
   Costin–Donninger–Glogić quasi-solution method), spectral gap ω0.
5. Riesz projections, linear stability on the stable subspace.
6. Nonlinear stability: Lipschitz bounds on the nonlinearity, Duhamel
   fixed point with Lyapunov–Perron correction / modulation of T (and
   boost parameters), statement of the theorem (data space, rate e^{-ωτ},
   ω < ω0), discussion of lower-order terms.
7. Numerical confirmation and honest list of gaps/caveats.
