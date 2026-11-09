# Sharp opening rates of single-strip subsolutions for the Muskat problem with viscosity jump

This folder holds a draft manuscript (`main.pdf`) on mixing solutions of the incompressible porous media (IPM) equation for two
fluids with different densities and different viscosities. It treats horizontally periodic interfaces and works within
*single-strip subsolutions*: one strip of constant intermediate density ρ_m around a moving central curve, opening linearly in time.

The paper builds on P. Luna-Velasco and F. Mengual, *Mixing zone dynamics for the Saffman–Taylor instability with viscosity jump*,
arXiv:2610.34504 ("LM26"). LM26 construct mixing solutions for fully unstable data (σ₀ < 0 everywhere) on asymptotically flat
interfaces. Credit for the ingredients:
- **LM26:** the Saffman–Taylor formulation (inverting in terms of σ rather than the vorticity, including the relation between γ
  and σ), the single-layer construction with its first- and second-order steps, and the flat single strip with arbitrary
  intermediate density.
- **Córdoba–Córdoba–Gancedo (Ann. of Math. 2011):** the periodic operator theory.
- **Mengual (Anal. PDE 2022):** the relaxation and the h-principle; Székelyhidi (2012) for A = 0.

## Results

Notation: A is the viscosity Atwood number, γ₀ is the vortex-sheet strength of the initial interface, and

    Φ(κ) = κ (2(1 − Aρ_m) κ + σ₀) − (1 − A²) γ₀² / (8 (1 + Aρ_m)).

- **Sufficient condition (Theorem 1.1).** Let z₀ be of class C^{6,γ}. Every opening κ ∈ C^{4,γ} with 0 < κ < κ_max (that is,
  κ > 0 and Φ(κ) < 0) is realised by a single strip, and the strip gives infinitely many mixing solutions.
- **Necessary condition (Theorem 1.2).** Let z₀ be of class C^{3,γ}. In every single-strip subsolution the centre of mass starts
  with the Muskat velocity and Φ(κ) ≤ 0. In particular κ = 0 on the pinch set {σ₀ ≥ 0, γ₀ = 0}.
- **The gap.** The two theorems are not a necessary-and-sufficient pair. Theorem 1.1 covers neither openings with Φ(κ) = 0
  somewhere nor openings that vanish at some points. The second case includes points where γ₀ ≠ 0, at which Φ(0) < 0 holds.
- **Pinch-free data (Corollary 1.3).** For z₀ of class C^{6,γ}, single strips with everywhere positive opening exist if and
  only if the pinch set is empty. For A ≠ 0 this includes partially Saffman–Taylor-stable data. The overturned example of
  Section 7 is pinch-free for A = 0 and, numerically, for A = 0.3, −0.6 and 0.9.

## What is new relative to LM26, and what is not

The audit in `NOVELTY.md` found that the revision of 7 Nov 2026 misdescribed LM26 and overstated its own novelty. The current text
states the novelty as follows.

New, but modest:
1. **Theorem 1.2.** Necessity within single-strip subsolutions of constant intermediate density. It is a short argument built
   from LM26's ingredients, and the flat case is in LM26 §3.2.
2. **An explicit vortex-sheet-dependent single-layer range for A ≠ 0.** It includes openings on pinch-free, partially
   Saffman–Taylor-stable data. It is latent in LM26's construction but not stated there.
3. **The optimal tilt δ\* for A ≠ 0 and general ρ_m.** It is sharp within single strips. Over LM26's own layer it gains a
   factor (2+A²)²/(4(1+A²)) ∈ [1, 9/8) on the γ₀² term, and nothing at A = 0.
4. **Curved strips with ρ_m ≠ A, and the horizontally periodic Hölder setting.** Both are routine.

Not new:
- **LM26's single layer is already tilted.** In the manuscript's variables it has D·τ₀ = −(1−A²)γ₀/(2+A²) at leading order,
  which equals δ\* when A = 0. With its exact hull condition evaluated in full, it satisfies the hull at t = 0 exactly when
  κ(2(1−A²)κ + σ₀) < (1−A²)γ₀²/(2(2+A²)²). LM26's stated range c̄ < 1/2 comes only from the last estimate of their §7.4, which
  bounds |Au + e₂| below by its normal component.
- **LM26's multilayers (§7.5) reach any c̄ < 1 for fully unstable data.** They open faster than single strips wherever |γ₀| is
  small.
- **Equal viscosities (Corollary 1.4).** For A = 0 the existence statements follow from LM26's construction. They are very
  likely covered by Förster–Székelyhidi (CMP 2018) for stable non-horizontal data and by Castro–Faraco–Mengual (Ann. PDE 2022)
  for turned interfaces. We have not compared our A = 0 statements with their precise admissible ranges.
- **Injection-driven viscous fingering (Corollary 1.5).** This is a change of frame and of time scale, and it applies equally
  to LM26's Theorem 1.1.
- **The pinch obstruction, qualitatively.** Mengual's abstract (Men22) already says that a pinch singularity in the
  relaxation prevents the two fluids from mixing wherever there is neither Rayleigh–Taylor nor vorticity at the interface.
  The pinch set {σ₀ ≥ 0, γ₀ = 0} is named after it; Remark 3.3(b) of the manuscript gives the precise relation.

Rigour caveats (see also `NOVELTY.md` §7):
- The time-regularity lemma and the C^{3,γ} regularity of G are only sketched (Appendix A.4).
- The h-principle (Men22; LM26 Thm 2.3) is applied to the 2π-periodic lift on ℝ², so the solutions live on ℝ² and need not be
  periodic; nothing is claimed on 𝕋×ℝ. (The 7 Nov revision asserted the h-principle on 𝕋×ℝ without proof.)
- Theorem 1.2 uses the uniform Plemelj–Privalov bounds of point (a) of Appendix A.4, which are not proved (Remark 6.3).
- Mixing is shown in space–time only; LM26 prove it in every time slice.
- The pinch-freeness of the example for A ≠ 0 is numerical, without validated error bounds.

## History

- `e8efc0f` (6 Nov 2026): first version. It was reconstructed **only from the one-paragraph abstract** of the talk
  *"Mixing solutions for unstable Muskat data with viscosity jump"* by Paula Luna-Velasco, joint work with Francisco Mengual
  (ICMAT, 12 November 2026), before arXiv:2610.34504 was read. That abstract already announced the key idea of inverting in terms
  of σ rather than the vorticity.
- `164466a` (6 Nov 2026): comparison with arXiv:2610.34504 (`COMPARISON.md`).
- `a1c40e3` (7 Nov 2026): revision written after reading LM26, titled *Vorticity-assisted mixing and sharp opening rates for the
  Muskat problem with viscosity jump*. It added §2.5 (the operator 𝒯, Prop. 2.10, γ₀ and the pinch set), the tilt, the
  γ₀-dependent range and Theorem 1.2. It also described LM26's second-order step wrongly.
- `66d5a96` (9 Nov 2026): novelty audit (`NOVELTY.md`). The current text corrects the title, abstract, credit note and body
  accordingly.

## Files

- `NOVELTY.md`: the audit of the 7 Nov revision against LM26. It covers the main finding, the uncredited overlaps, the
  provenance, the outside literature, per-claim verdicts, what remains new, and rigour caveats.
- `COMPARISON.md`: detailed comparison with LM26, including the original comparison of the abstract-only reconstruction with
  the paper.
- `main.pdf`: compiled manuscript.
- `main.tex`, `sec-*.tex`, `app-cauchy.tex`, `biblio.tex`, `fig-turned.pdf`: LaTeX sources. Build with `latexmk -pdf main.tex`.
- `checks/`: Python scripts (numpy, sympy, matplotlib).
  - `check_sigma.py`: on a test graph, the Saffman–Taylor and vorticity formulations give the same normal velocity
    (difference ~1e-14). It also checks ∫σ ds = −4π and σ = −2(n₂ + AU).
  - `check_algebra.py`: symbolic checks of the first-order limits, a·n and the t = 0 hull identity, the Rankine–Hugoniot
    conditions, the mass identity, the flat example and the laminate lemma.
  - `check_vorticity.py`: the identity γ = −2(τ₂ − A𝒯σ)/(1−A²) (LM26 (3.7b) in our normalisation), and the stable fraction,
    min |γ₀| and min κ_max for the overturned example.
  - `check_sharp.py`: the sharp-range algebra: the completed square in the t = 0 hull identity, the optimality of δ\*, and κ_max
    as the root of Φ.
  - `check_hull_equiv.py`: the relaxed set 𝒰_A equals Mengual's set O (written (O) in LM26), with
    D·(D−2a) = 4(1−A²) g / ((1−ρ²)²(1+Aρ)), where g is LM26 (A.1).
  - `compare_hulls.py`: the same equivalence, as sign agreement on 200,000 random states.
  - `check_lm26_layer.py`: LM26's own single layer in the manuscript's variables (`NOVELTY.md` §1). (a) LM26 (4.14) is the
    disc |D−a| < |a|; (b) its tilt δ_LM = [2(1+A²)/(2+A²)] δ\*, its exact-hull range at t = 0 and the gain factor in [1, 9/8);
    (c) min κ_max against LM26's exact-hull bound on the overturned example; (d) M = −2t∂_α(κδ_LM) + O(t²) on LM26's layer.
  - `make_figure.py`: generates `fig-turned.pdf` (from this folder: `python3 checks/make_figure.py fig-turned.pdf preview.png`).
