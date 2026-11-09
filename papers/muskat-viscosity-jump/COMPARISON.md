# Reconstruction vs. the actual paper

## Correction (9 Nov 2026)

The 7 Nov table below, and the paragraph under it on how the gain arises, misdescribed LM26 and overstated what the
revision adds. The verified account is the audit in `NOVELTY.md` (9 Nov 2026). Where this file and `NOVELTY.md` differ,
`NOVELTY.md` holds. In short:
- **LM26's single layer is already tilted.** Their second-order step cancels the `O(t)` part of `∂_αG₊ − ∂_αG₋`, a
  condition on their potential `g`, not on our mismatch `M`. In our variables it gives the tilt
  `δ_LM = −(1−A²)γ₀/(2+A²) = [2(1+A²)/(2+A²)]δ*`, which is the optimal tilt `δ*` when `A = 0`.
- **Their hull condition already contains the vortex sheet.** (4.14), `2|∇g| < |Au+e₂|`, is our disc `|D−a| < |a|`. Their
  stated range `c̄ < 1/2` comes only from the lower bound `|Au+e₂| ≥ (Au+e₂)·τ^⊥ = −σ₀/2` in their §7.4. With the full
  modulus and `c` free, their unchanged layer satisfies the strict hull at `t = 0` iff
  `κ(2(1−A²)κ + σ₀) < (1−A²)γ₀²/(2(2+A²)²)`, so it also opens on pinch-free, partially Saffman–Taylor-stable data. This is
  latent in LM26 (a few-line change of their §7.4), not stated there.
- **Uncredited overlaps.** Prop. 2.10 is LM26 (3.7b) with their §3.1.3, and the `A = 0` formulas of Example 2.12
  (`σ₀ = −2n₀,₂`, `γ₀ = −2τ₀,₂`, fully unstable ⟺ graph) are LM26 §§1.2 and 3.1.3. In Lemma 5.2, `a^{(0)}·n₀` is LM26
  (7.11) with the computation of Prop. 7.4, and `a^{(0)}·τ₀` follows from (3.7b). Lemmas 5.5–5.6 are LM26's second-order
  step (Lemma 5.3, (5.21), (8.41)) with a nonzero target, solvable by LM26 Prop. 6.1 (periodic case: CCG11). `𝒰_A` is the
  set `O` of Mengual (Men22), written (O) in LM26, and `g` is LM26 (A.1) (Men22, Lemmas 4.1–4.2). The flat strip with free
  intermediate density is LM26 §3.2; the flat continuous profile (Example 3.6) and the untilted multilayer limit (Sec. 8.1)
  are Men22 and LM26 §§3.2 and 7.5. The periodic operator theory is that of Córdoba, Córdoba and Gancedo (CCG11). The
  injection corollary applies equally to LM26 Thm 1.1.
- **What remains new** is listed in `NOVELTY.md` §6. The optimal tilt `δ*` gains a factor `(2+A²)²/(4(1+A²)) ∈ [1, 9/8)` on
  the `γ₀²` term over LM26's own layer, and nothing at `A = 0`. LM26's multilayers with several intermediate strips open
  faster than single strips for fully unstable data wherever `|γ₀|` is small.

The 7 Nov table has been corrected in place. The original comparison of 6 Nov is kept below; its wrong statements are
marked *[Corrected 9 Nov 2026: …]* and later additions *[9 Nov: …]*. Lemma numbers refer to the revision of 9 Nov.
`checks/check_lm26_layer.py` verifies the formulas quoted for LM26's layer.

## Revision of 7 Nov 2026 (table corrected 9 Nov 2026)

The manuscript has been revised to go beyond Luna-Velasco & Mengual (arXiv:2610.34504, "LM26") and is now titled
*Vorticity-assisted mixing and sharp opening rates for the Muskat problem with viscosity jump*. LM26 is cited as the source of the
fully unstable result, the Saffman–Taylor formulation and the single-layer construction. *[Corrected 9 Nov 2026: the gain over
LM26 is much smaller than "go beyond" suggested; see the Correction above. The revision of 9 Nov is titled* Sharp opening rates
of single-strip subsolutions for the Muskat problem with viscosity jump*.]*

| | LM26 | Manuscript (as corrected 9 Nov) |
|---|---|---|
| Admissible data | Stated: fully unstable, `σ₀ < 0` everywhere (Thm 1.1). With the exact hull (4.14), their construction also opens where `σ₀ ≥ 0` and `γ₀ ≠ 0`: latent, not stated (a few-line change of §7.4, with `c − c_∞ ∈ H^{k−1}` and `inf c > 0`) | **pinch-free**: no point with `σ₀ ≥ 0` and `γ₀ = 0` (Cor. 1.3); includes partially Saffman–Taylor *stable* interfaces. The existence half is latent in LM26's layer; the "only if" half is Thm 1.2 |
| Single-layer opening range | Stated: `c̄ < 1/2`, i.e. `2(1−A²)κ < −σ₀` ((7.13)), from the lower bound `\|Au+e₂\| ≥ (Au+e₂)·τ^⊥ = −σ₀/2` in §7.4; `c` is free in §§5–7.3 ((5.1)). Exact hull of the unchanged layer at `t = 0`: `κ(2(1−A²)κ + σ₀) < (1−A²)γ₀²/(2(2+A²)²)` (not stated) | `κ(2(1−Aρ_m)κ + σ₀) < (1−A²)γ₀²/(8(1+Aρ_m))` (Thm 1.1), sharp within single strips (Thm 1.2). For `ρ_m = A` the `γ₀²` coefficient is that of LM26's layer times `(2+A²)²/(4(1+A²)) ∈ [1, 9/8)`, equal at `A = 0` |
| Multilayers | `2N` interfaces (`2N−1` intermediate strips with the densities `θ_j` of (7.14); the middle one has `θ_N = A`, and for `N = 1` the levels are −1, A, 1, a single strip with `c̄ < 1/2`) reach any `c̄ < 1` for fully unstable data (§7.5); this exceeds the single-strip `κ_max` wherever `\|γ₀\|` is small | single strips only |
| Sharpness | Curved data: not addressed; LM26 say their main goal is "to construct mixing solutions rather than optimize the opening speed" (§7.5). Flat data: §3.2 finds the admissible openings of its single-strip ansatz (hull `\|γ\| < 1/2`, mass balance) and optimises `θ_mix` | **necessary** within single-strip subsolutions of constant intermediate density (Thm 1.2): a short argument built from LM26's ingredients, whose flat case is LM26 §3.2. For LM26's own layer the conclusions follow from their formulas; the transfer of Thm 1.2 to an asymptotically flat class is expected, not proved (Steps 1 and 4 use periodic tools) |
| Second-order step | `z^{(2)}` cancels the `O(t)` part of `∂_αG₊ − ∂_αG₋` (§3.3.4, Lemma 5.3), so `∂_λG = O(t²)` and `∇g` is tangential at leading order (Props. 7.2–7.3, (7.9)). In our variables: tilt `δ_LM = −(1−A²)γ₀/(2+A²) = [2(1+A²)/(2+A²)]δ*`, equal to `δ*` at `A = 0`, and `M = −2t∂_α(κδ_LM) + O(t²)` | the same step (affine in the acceleration through `I − 2A𝒦₀ ↔ I − AP`, LM26 Lemma 5.3, (5.21), (8.41); solvable by LM26 Prop. 6.1, periodic case CCG11), with its target chosen to give the optimal tilt `δ* = −(1−A²)γ₀/(2μ_m)` (Lemmas 5.5–5.6): an optimisation of LM26's normalisation, not a new mechanism. The untilted strip `δ = 0` is a different subsolution, not LM26's |
| Equal viscosities, A = 0 | Not stated. Since `\|Au+e₂\| = 1`, (4.14) and (7.9) give the range `0 < κ < (1+n₀,₂)/2` for their layer directly | mixing along the whole of any overturned interface whose tangent never equals −e₁ (Cor. 1.4), which also follows from LM26's construction. Very likely covered by FS18 (stable non-horizontal data) and CFM22 (turned interfaces); we have not compared with their precise admissible ranges |
| Injection-driven viscous fingering (no gravity) | Not stated; follows from LM26 Thm 1.1 by the same change of frame (fully unstable data) | Cor. 1.5, via the effective gravity `G + AV`: a change of frame and of time scale that applies equally to LM26 Thm 1.1; new only through the pinch-free range |
| Intermediate density | Flat strip: `θ_mix` free (§3.2), and `θ_mix = A` is chosen because it maximises `h₊h₋`. Curved single layer: `θ = A`. Multilayers: the levels `θ_j` of (7.14), which include `θ_N = A` and, for `N ≥ 2`, other values; (4.6) and (4.10) are written for general `θ` and reused with `θ = θ_j` ((7.17)) | curved single strips with any `ρ_m ∈ (−1,1)` (routine) |
| Vortex strength vs σ | (3.7b), left relation: `(1−A²)ϖ = ⟦θ⟧(AQ − i)·∂_αz`; recovering `ϖ` from `σ`, and the operators `I ± AP` (§3.1.3) | `γ = −2(τ₂ − A𝒯σ)/(1−A²)` (Prop. 2.10): the same identity in our normalisation (`ϖ = −\|∂_αz\|γ`, `⟦θ⟧ = −2`, `Q·τ ↔ 𝒯σ`), written after reading LM26 |
| Relaxed set | `O` of Men22, written (O) (§2); `Ō = K^Λ = R` (Thm 2.4; `K^{lc,Λ} = Ō` is due to Men22) | `𝒰_A = O`, a reformulation at fixed Darcy field: `D·(D−2a) = 4(1−A²) g / ((1−ρ²)²(1+Aρ))`, where `g` is the function of LM26 (A.1) (Men22, Lemmas 4.1–4.2) |
| Mixing property | in every time slice (Def. 2.2) | in space–time only |
| Setting | asymptotically flat, `z° − (α,0) ∈ H⁵`, chord–arc | horizontally periodic, `C^{6,γ}`; the periodic operator theory is CCG11's, as LM26 note (§3.1.3) |

*History of this table.* The 7 Nov version stated, in the LM26 column: a range "(σ only)"; "`z^{(2)}` chosen to cancel the
mismatch (no tilt)"; "not covered" for injection; "—" for the γ–σ identity and for `A = 0`; "`θ_mix = A`" for the intermediate
density; and "not addressed" for sharpness. It listed "`𝒰_A`; proved equal to `O`" and the free `ρ_m` among the results
beyond LM26, and it said that Thm 1.2 also applies to LM26's single layer. These entries were wrong or incomplete; see
`NOVELTY.md` §§1, 2 and 5.

**How the γ₀ term arises (corrected 9 Nov 2026).** The relaxed flux `D` must lie in the disc `|D−a|<|a|`. Its normal part is
fixed by the Rankine–Hugoniot conditions. Its tangential part (the tilt) is set by the `O(t)` term of the flux mismatch, which
the second-order acceleration controls. At `t = 0` the centre `a` has normal component `−σ₀/2 − A²(1−ρ_m²)κ/μ_m` and
tangential component `−(1−A²)γ₀/(2μ_m) = δ*`, so the admissible range is governed by `|a|`, which combines the
Saffman–Taylor function σ₀ and the vortex sheet γ₀ (Lemma 5.2; for `ρ_m = A` the normal component is LM26 (7.11) with the
computation of Prop. 7.4, and the tangential one follows from (3.7b)). This account of the two components of `D` follows
LM26's own description of their construction (§§3.3.1 and 3.3.4, (4.11)–(4.12)), with their potential `g` in place of `Ψ`.
*[Corrected 9 Nov 2026: the 7 Nov text said here that "LM26
make that term vanish; we instead tilt `D` toward the centre of the disc". That was wrong. LM26's normalisation gives the tilt
`δ_LM`, towards the centre and exactly onto it at `A = 0`. Their hull condition (4.14) already uses the full modulus
`|Au+e₂|`, which carries γ₀ through `(Au+e₂)·τ₀ = δ_LM` at leading order; γ₀ is lost only in their final lower bound (§7.4).
The manuscript adds (i) the evaluation of the hull exactly instead of through its normal component and (ii) the optimal
choice `δ*` of the target of LM26's second-order condition, which can be prescribed freely because `I − AP` is invertible.
Step (ii) multiplies the γ₀² term by a factor in `[1, 9/8)`, and by 1 at `A = 0`.]*

**Checks (in `checks/`):**
- `check_vorticity.py`: verifies the identity for γ, i.e. LM26 (3.7b) in our normalisation (a verification, not a new result),
  with maximal absolute residual ≤ 2·10⁻¹³. On the overturned interface `(α−1.5 sin α, 0.8 sin α)` and for
  `A = 0, 0.3, −0.6, 0.9`, `σ₀ ≥ 0` on about 27%, 26%, 24%, 16% of the parameter interval (15%, 15%, 14%, 9% of the
  arclength; the arclength fractions are printed by `check_lm26_layer.py`), with `min |γ₀| ≈ 1.70, 1.74, 1.88, 2.33` there.
  (Both scripts count the nodes of the 512-point grid, so the fractions have a resolution of about 0.2%.) So the data are pinch-free (for `A ≠ 0` numerically, without
  validated error bounds), and `κ_max ≥ 0.22`. *[Corrected 9 Nov 2026: the 7 Nov text said "16–27% of the curve", "|γ₀| ≥ 1.7"
  and "the old range is empty on that part". The percentages are fractions of the parameter interval, `min |γ₀|` is 1.696 at
  `A = 0`, and only the range stated by LM26 is empty there: their own layer with the exact hull opens there as well.]*
- `check_lm26_layer.py` (added 9 Nov): LM26's single layer in our variables. It checks `D = (Au+e₂) − 2∇^⊥g`,
  `a = (e₂ − A∇p)/(1+A²)` and `D·(D−2a) = (4|∇g|² − |Au+e₂|²)/(1+A²)`; the tilt `δ_LM`; the exact range of their layer; the
  identity of Lemma 5.2; and the relation between `M` and `∂_αG₊ − ∂_αG₋`. On the overturned example the minima over the
  interface are `κ_max = 0.2350, 0.2297, 0.2265, 0.2959` against `0.2350, 0.2294, 0.2224, 0.2778` for LM26's layer with the
  exact hull (`A = 0, 0.3, −0.6, 0.9`).
- `check_sharp.py`: the sharp-range algebra.
- `check_hull_equiv.py`: the identity `D·(D−2a) = 4(1−A²) g / ((1−ρ²)²(1+Aρ))` with a positive factor, so `𝒰_A = O`, the set
  of Men22 (a reformulation, not a new result). `compare_hulls.py` tests the same identity on 200,000 random states.
- `make_figure.py`: produces Figure 1.

**Caveats.**
- For `A = 0` our statements are very likely covered by earlier work: see Förster–Székelyhidi (CMP 2018) for stable
  non-horizontal data and Castro–Faraco–Mengual (Ann. PDE 2022) for turned interfaces. We have not compared our `A = 0`
  statements with their precise admissible ranges; the full texts were not reachable, and `NOVELTY.md` §4 records what the
  abstracts and search snippets say. *[Corrected 9 Nov 2026: the 7 Nov text described CFM22 only as constructing mixing zones
  localised near the unstable part of turned interfaces. According to its abstract and search snippets, it opens the mixing
  zone slightly inside the stable regime, and it reportedly identifies "fully stable points" (horizontal and stable) as the
  obstruction.]*
- For `A ≠ 0`, Mengual's abstract (Men22) already says that a pinch singularity in the relaxation prevents mixing wherever
  there is neither Rayleigh–Taylor nor vorticity at the interface.
- Rigour: the time-regularity lemma (Lemma 5.4), the `C^{3,γ}` regularity of `G` and the differentiability at `t = 0` of the
  coefficient maps in Lemma 5.5 are only sketched (Appendix A.4), as in the first version. LM26 prove the asymptotically flat
  analogues (Lemmas 8.4–8.6, Prop. 7.1; for the second time derivative only the weaker bound on spatial primitives of
  Lemma 8.6 and Prop. 7.2, not a Lipschitz bound on `∂_tF^±`), so the asymptotically flat analogue of the existence half rests
  on their estimates.
  The h-principle is applied to the periodic lift on `ℝ²`. Mixing is shown in space–time only.

---

## Original comparison (reconstruction from the abstract vs. LM26)

**Reconstruction:** `main.pdf` in this folder (20 pp.). It was written from the seminar abstract only and committed
(`e8efc0f`) before the actual paper was opened. *[9 Nov: this is the `main.pdf` of commit `e8efc0f`; the `main.pdf` now in
this folder belongs to a later revision. The seminar abstract already announced the key idea of inverting in terms of σ
rather than the vorticity, so the independence concerns details, not the idea.]*

**Actual paper:** P. Luna-Velasco and F. Mengual, *Mixing zone dynamics for the Saffman–Taylor instability with viscosity jump*,
arXiv:2610.34504v1 (25 Oct 2026), 45 pp.

## 1. Verdict

The reconstruction recovers the paper's architecture almost one-to-one:

- the same reformulation (pressure / Darcy field instead of the vorticity);
- the same new Muskat formulation;
- the same type of subsolution (one intermediate density, opening along the initial normal);
- the same second-order expansion in time with the same two compatibility conditions; *[Corrected 9 Nov 2026: the first-order
  conditions agree, the second-order ones do not. LM26 cancel the `O(t)` part of `∂_αG₊ − ∂_αG₋`, the reconstruction the
  `O(t)` part of `M`. On LM26's layer `M = −2t∂_α(κδ_LM) + O(t²)`, so the two order-`t` conditions differ by the
  acceleration-independent term `−2∂_α(κδ_LM)`.]*
- the same Fredholm/energy inversion of `I − A·P`; *[9 Nov: LM26 credit the argument to CCG11, Sec. 4.]*
- the same relaxed set (shown numerically to coincide). *[9 Nov: the coincidence is an exact identity (`check_hull_equiv.py`),
  and the set is the set `O` of Men22, written (O) in LM26.]*

With the intermediate density chosen as `ρ_m = A`, the reconstruction gives exactly the paper's single-layer mixing zone and its admissible opening range. *[Corrected 9 Nov 2026: only the `O(t)` geometry of the strip and the range stated by LM26
coincide. The subsolutions differ: LM26's layer has the tilt `δ_LM`, the reconstruction's strip `δ = 0`. Their central curves
differ at `O(t²)` (different accelerations) and in the tangential part of `z^{(1)}`. And LM26's construction with the exact
hull admits the larger, γ₀-dependent range `κ(2(1−A²)κ+σ₀) < (1−A²)γ₀²/(2(2+A²)²)`.]*

The paper goes further in four places:
- **Exact relaxation:** it proves the relaxation is exact. *[9 Nov: LM26 prove `K^Λ = R = Ō` (Thm 2.4); the identity
  `K^{lc,Λ} = Ō` is due to Men22, as LM26 state.]*
- **Multilayers:** it pushes the opening rate to the optimal range. *[Corrected 9 Nov 2026: to any `c̄ < 1` (§7.5). "Optimal"
  is not LM26's claim: they say their main goal is to construct mixing solutions rather than optimise the opening speed.]*
- **Setting:** it works in the asymptotically flat, `H^5` setting. *[9 Nov: the periodic operator theory used by the
  reconstruction is due to CCG11, to which LM26 attribute it.]*
- **Proofs:** it gives fully explicit second-order formulas and about 12 pages of operator estimates.

On those steps the reconstruction is shorter and partly at sketch level.

## 2. Side-by-side

| Item | Reconstruction | Luna-Velasco–Mengual |
|---|---|---|
| Title | Mixing solutions for unstable Muskat data with viscosity jump (talk title) | *Mixing zone dynamics for the Saffman–Taylor instability with viscosity jump* |
| Normalisation | `ρ=±1`, `μ=1+Aρ`, `μu=−∇p−ρe₂` | `θ=±1`, `u+Aθu=−∇p−θe₂` (identical) |
| Setting | Horizontally periodic (`T×R`) *[9 Nov: the estimates for `P` and the bounds for `(I−ξP)⁻¹`, `\|ξ\| ≤ 1`, in this setting are due to CCG11, as LM26 note in §3.1.3]* | Asymptotically flat on `R` (periodic only in their Figure 1) |
| Regularity | `z₀ ∈ C^{6,γ}` (Hölder) | `z₀−(α,0) ∈ H⁵`, chord–arc |
| Saffman–Taylor function | `σ = [∂ₙp]`, per arclength | `σ = −⟦∇p⟧·∂_αz^⊥`, per unit `α` (same object, `σ_theirs = \|∂_αz\| σ_mine`) |
| Muskat in σ-form | `(I−2A𝒦)σ = −2n₂`, `U = −𝒦σ` | `(I−AP)σ = −2∂_αz₁`, `(∂ₜz+Q)·∂_αz^⊥ = 0`; "(Mσ) seems to be new" |
| Dual (vorticity) form | `(I+2A𝒦)γ = −2τ₂` (Remark 2.9) | `(I+AP)ϖ = −⟦θ⟧∂_αz₂`: "both problems reduce to inverting operators of the form I ± AP" |
| Why σ helps | `v = μu` has explicit curl/div; vorticity of `u` involves `curl m` *[Corrected 9 Nov 2026: the Darcy field is `v = u + Am = −∇p − ρe₂`, which equals `μu` only on `K`]* | "pressure remains harmonic in each region, whereas the mixing zone generally carries interior vorticity", `(1+A²/2) curl u = AΔg` |
| Relaxed set | `{\|ρ\|<1, \|D−a\|<\|a\|}`; shown ⊂ second-order laminates *[9 Nov: this is the set `O` of Men22 sliced at fixed Darcy field, and the inclusion re-derives part of Men22's result]* | `O = {\|θ\| < 1, g(z) < 0}`, with `g(z) := (1−Aθ)\|m−θu\|² + (1−θ²)(m−θu)·(Au+e₂)` (LM26 (A.1) characterises `Ō` by `\|θ\| ≤ 1, g ≤ 0`); proves `K^{lc,Λ} = K^Λ = R = Ō` (Appendix A, with AI-suggested Λ-affine functions) *[9 Nov: LM26 credit `K^{lc,Λ} = Ō` and the characterisation by `g` (A.1) to Men22, Lemmas 4.1–4.2; their new part is `K^Λ = R = Ō` (Thm 2.4)]* |
| h-principle | Imported from Mengual (APDE 2022) | Same source, their Theorem 2.3 |
| Subsolution | One strip, density `ρ_m ∈ (−1,1)` arbitrary; speeds `(1±ρ_m)κ` | One strip, `θ_mix = A`; speeds `c/(1∓A)` *[9 Nov: for curved data. The flat strip of §3.2 keeps `θ_mix` free and chooses `θ_mix = A` because it maximises `h₊h₋`; the multilayers use the levels `θ_j` of (7.14), which include `θ_N = A` and, for `N ≥ 2`, other values]* |
| Opening rate | `0 < κ < −σ₀/(2(1−Aρ_m))` | `c = −c̄ σ°/\|∂_αz°\|`, `c̄ < 1/2` (single layer) *[9 Nov: this is the stated range, which comes from the lower bound `\|Au+e₂\| ≥ (Au+e₂)·τ^⊥` in §7.4. `c` is free in §§5–7.3, and with the exact hull their layer admits `κ(2(1−A²)κ+σ₀) < (1−A²)γ₀²/(2(2+A²)²)`]* |
| Strip momentum | Stream function `Ψ`, `D = ∇^⊥Ψ`, linear in λ between boundary values from Rankine–Hugoniot | Potential `g`, `m = θu − (1−θ²)/(1−Aθ)(∇^⊥g + ½(Au+i))`, linear in λ between boundary values from Rankine–Hugoniot *[9 Nov: in our variables `D = (Au+e₂) − 2∇^⊥g` exactly, so the boundary fluxes are `F^± = (AU^±+e₂)·∂_αz_±^⊥ − 2∂_αG_±`]* |
| Central curve | `Z = z₀ + (tU₀ + ½t²a)n₀ + h(t)e₂` | `z = z° + t z^{(1)} + ½t² z^{(2)}` |
| Order `t⁰` | Flux mismatch `M(0)=0` ⇔ Muskat velocity | `∂_αG₊ = ∂_αG₋` at `t=0` ⇐ `z^{(1)} = −Q₀[σ^{(0)}]` (LM26 Lemma 5.2; only the normal component of `z^{(1)}` enters, and LM26 do not state the converse) *[9 Nov: the same condition; on LM26's layer `M(·,0) = −2(∂_αG₊ − ∂_αG₋)(·,0)`]* |
| Order `t¹` | `(I−2A𝒦₀)a = R₁` (remainder only shown to exist and be `C^{3,γ}`) | Explicit `(I−AP)σ^{(1)} = …`, `z^{(2)}` from `(5.16)` with explicit remainders `r_σ, r_G, R_a` *[Corrected 9 Nov 2026: the two targets differ. LM26 impose `∂_t(∂_αG₊ − ∂_αG₋)(0) = 0` ((5.19)), which in our variables reads `∂_tM(0) = −2∂_α(κδ_LM)`; the reconstruction imposed `∂_tM(0) = 0`]* |
| Key smallness | `Ψ⁺−Ψ⁻ = ∫M = O(t²)` | `∂_λG = ½(G₊−G₋) = O(t²)` *[Corrected 9 Nov 2026: the two are not equivalent. LM26's condition makes the normal component of `∇g` of order `t`. On their layer, in our variables, `Ψ⁺ − Ψ⁻ = −tκ∫₋₁¹(Au+e₂)·τ₀ dλ − 2(G₊ − G₋) = −2tκδ_LM + O(t²)`, and `M = −2t∂_α(κ(Au+e₂)·τ₀) − 2(∂_αG₊ − ∂_αG₋) + O(t²)`. The correspondence `M ↔ ∂_αG₊ − ∂_αG₋` holds only at order `t⁰`]* |
| Mass / primitive | Mass identity `∫_T M = −2Q/(1−ρ_m²)` plus a vertical shift `h(t) = O(t³)` | Primitive estimates on `R` (Lemma 8.6 *[9 Nov: with Prop. 7.2 and Lemma 8.2]*), no shift needed |
| Inversion of `I − A·P` | Riesz–Schauder + harmonic-energy argument (Hölder) | Riesz–Schauder + harmonic-energy argument on the adjoint (`L²`, then `H^j`) *[9 Nov: "the classical harmonic-energy argument; cf. [15, Section 4]", i.e. CCG11, who are also credited for the smoothing of `P` (Lemmas 3.1 and 5.2 there)]* |
| Multilayer / optimal rate | Formal remark (§6.2 *[9 Nov: of `e8efc0f`; now §8.1]*): speeds up to `−σ/(1−A)`, `−σ/(1+A)` | Proved (§7.5): `2N` layers reach any `c̄ < 1`, i.e. exactly those speeds *[Corrected 9 Nov 2026: `2N` interfaces (`2N−1` intermediate strips) reach any `c̄ < 1`; the speeds `−σ/(1∓A)` are the supremum (`c̄ → 1`), approached but not reached. LM26 do not call them optimal]* |
| Flat case | Continuous profile with flux `θ(1−ρ²)/(1−Aρ)`, fronts `2θ/(1∓A)` | One strip: `h_± = 1/(1∓A)`; continuous profile flux `F_A = (1−θ²)/(1−Aθ)` cited from Mengual 2022 *[9 Nov: §3.2 obtains `h_± = 1/(1∓A)` from a flat strip with `θ_mix` free (hull `\|γ\| < 1/2`, mass balance). For the flat interface, LM26 Thm 1.1 with the multilayers of §7.5 and `c̄ = θ` gives mixing zones equal to those of the reconstruction's profile, and LM26 identify the `N → ∞` limit of their profiles with the Otto profile with time scaled by `c̄`, i.e. this profile]* |
| Extras | Numerical/symbolic checks (`checks/`) | Vortex-blob simulation (Fig. 1); remark on spontaneous stochasticity |

## 3. What matched

1. **The key idea.** The paper says it inverts "in terms of the Saffman–Taylor function rather than the vorticity" because the pressure (equivalently my `v = μu`) is harmonic in each region. That is exactly the reconstruction's rationale. Both papers also notice the duality of the two formulations (`I ± A·P`). *[Corrected 9 Nov 2026: the Darcy field of a
   mixing-zone state is `v = u + Am = −∇p − ρe₂`, which equals `μu` only on `K`. The idea was announced in the seminar
   abstract from which the reconstruction was written, so this agreement is not independent. The rationale and the duality
   are LM26's (Secs. 3.1.3–3.1.4 and 3.3.5 of their preprint of 25 Oct 2026), and LM26 attribute the periodic estimates for
   `P` and the inverse bounds for `I − ξP` to CCG11.]*
2. **Same equations.** The σ-equations agree once normalisations are matched (`P ↔ 2𝒦`, `σ_theirs = |∂_αz| σ_mine`). Their boundary densities `σ_a^{(0)} = (1−aA)σ^{(0)}/2 − aAc|∂_αz°|` are the negatives of the source parts of my densities `f^±` at `t=0` with `ρ_m = A`.
3. **Same relaxation.** On 200,000 random states, my condition `D·(D−2a) < 0` and their `g < 0` agree every time; the ratio of the two quantities depends only on `(θ, A)`. Inside the strip with `θ = A` the equivalence `|D−a|<|a| ⇔ 2|∇g| < |Au+e₂|` is exact algebra. *[9 Nov: the ratio is `4(1−A²)/((1−θ²)²(1+Aθ)) > 0`, now asserted in
   `check_hull_equiv.py` and `compare_hulls.py`; the common set is the set `O` of Men22.]*
4. **Same geometry and admissible range.** With `ρ_m = A` and `κ = c/(1−A²)`, my curves `z^± = Z ± t(1±A)κ n₀` are their `z_± = z ± t c τ^⊥/(1∓A)`. My range `(1−A²)κ < −σ₀/2` is their `c̄ < 1/2`. The positivity margins vanish at the same endpoints. *[Corrected 9 Nov 2026: the first-order geometry coincides, and `c̄ < 1/2` is the range LM26
   state. But their exact margin `|Au+e₂| − 2|∇g|` does not vanish at `c̄ = 1/2` when `γ₀ ≠ 0`; only their lower bound
   `|Au+e₂| ≥ (Au+e₂)·τ^⊥` (§7.4) does. The reconstruction's untilted strip is not LM26's subsolution.]*
5. **Same proof skeleton.** Both papers:
   - follow the time expansion of Förster–Székelyhidi;
   - interpolate a scalar potential linearly across the strip;
   - match first order with the Muskat velocity;
   - match second order by solving the same operator equation again;
   - reduce the tangential part of the relaxed flux to a primitive of size `O(t²)`. *[Corrected 9 Nov 2026: LM26 make
     `∂_λG = O(t²)`, i.e. the normal component of `∇g` of order `t` (Prop. 7.2). In our variables this fixes the tangential
     part of `D` at `δ_LM + O(t)`, not at `O(t)`. The second-order targets differ accordingly; see the rows "Order `t¹`" and
     "Key smallness" above.]*
6. **Optimal rates.** My formal "continuous profile" remark predicted fronts at `−σ/(1−A)` and `σ/(1+A)` relative to the central curve. These are exactly the limits of their multilayer construction ("converge to the viscosity-dependent Otto profile"). *[Corrected 9 Nov 2026: "predicted" carries no priority. The remark was written on 6 Nov
   (`e8efc0f`), after LM26 v1 (25 Oct), and only conjectured what LM26 §7.5 proves. The quoted "converge to the
   viscosity-dependent Otto profile" is LM26's `N → ∞` limit at fixed `c̄`, whose fronts are `c̄` times these speeds. The
   fronts above are the supremum `c̄ → 1`, approached but not reached, and LM26 do not call these rates optimal.]*

## 4. Where they differ

- **Setting and regularity.**
  - The reconstruction is periodic and Hölder (`C^{6,γ}`); the paper is asymptotically flat and Sobolev (`H⁵`).
  - The periodic choice forced an extra device: the mass identity and the vertical shift `h(t)`, to keep the strip's stream function periodic. *[Corrected 9 Nov 2026: to keep `D = ∇^⊥Ψ` periodic; `Ψ` must have equal increments
    `∫_T F⁺ = ∫_T F⁻` over a period on both boundaries.]*
  - On `R` the paper controls primitives directly *[9 Nov: Prop. 7.2, Lemmas 8.2 and 8.6]*.
- **Intermediate density.**
  - The reconstruction allows any `ρ_m`. This is a one-parameter family of asymmetries, e.g. `ρ_m → 1` gives the optimal upper speed but a vanishing lower one. *[Corrected 9 Nov 2026: stale and unproved. This used the 6 Nov range
    `κ < −σ₀/(2(1−Aρ_m))`; `ρ_m = ±1` is excluded, no optimality in `ρ_m` is proved, and for the flat interface the statement
    follows from the mass balance of LM26 §3.2.]*
  - The paper fixes `θ_mix = A`, which maximises the flat-case product `h₊h₋`. *[9 Nov: LM26 §3.2 derive this choice for a
    flat strip with `θ_mix` free. Their curved single layer uses `θ = A`, and their multilayers use the levels `θ_j` of (7.14),
    which include `θ_N = A` and, for `N ≥ 2`, other values.]*
- **Explicitness of the second-order step.**
  - The paper computes everything explicitly: `σ^{(1)}`, `z^{(2)}`, and the remainders `R_a, r_a, r_σ, r_G`.
  - The reconstruction (Lemma 5.4 *[Corrected 9 Nov 2026: Lemma 5.5]*) only proves that `∂ₜM(0)` is affine in the acceleration, with linear part `−(4|∂_αz₀|/(1−ρ_m²))(I−2A𝒦₀)`. It argues the remainder exists and is `C^{3,γ}` without writing it. *[9 Nov: this is
    LM26's second-order step (Lemma 5.3, (5.21), (8.41)), solvable by LM26 Prop. 6.1 (periodic case: CCG11).]*
- **Relaxation.**
  - The reconstruction only shows its set lies in the lamination hull; this is enough for the h-principle. *[9 Nov: this
    re-derives one inclusion of Men22's identity `K^{lc,Λ} = Ō`.]*
  - The paper proves the full identity `K^{lc,Λ} = K^Λ = R = Ō`, including the bounded version
    `K_M^{lc,Λ} = K_M^Λ = R_M = Ō_M` (Thm 2.5). *[9 Nov: LM26 credit
    `K^{lc,Λ} = Ō` and its bounded version to Men22; their new part is `K^Λ = R` (Thms 2.4–2.5).]*
- **Opening rates.**
  - The reconstruction proves only the single-strip range.
  - The paper proves `c̄ < 1` via finitely many layers (§7.5, presented briefly, with details in Lemma 8.8).

## 5. Rigor of the reconstruction, compared honestly

- **At the paper's level:**
  - the σ-formulation and invertibility (Prop. 2.6–2.7; checked numerically against the vorticity formulation to `~1e−14`) *[9 Nov: the periodic operator theory behind them is CCG11's]*;
  - the laminate computation and the relaxed set;
  - the algebra of the ansatz: jump conditions, mass identity, `t=0` limits, hull at `t=0` (all checked symbolically in `checks/`).
- **Thinner:**
  - the time-regularity lemma (5.3) and the regularity of the remainder `G` (Lemma 5.4) are a structured sketch in the appendix (derivative identities plus Plemelj–Privalov), not full estimates; *[Corrected 9 Nov 2026: Lemmas 5.4 and
    5.5. The differentiability at `t = 0` of the coefficient maps in Lemma 5.5 is also only sketched (Appendix A.4).]*
  - the tubular-neighbourhood lemma (4.1) is a sketch *[9 Nov: its counterparts are LM26 Lemma 4.1 and CFM22 Lemma B.1]*.
  
  The paper spends §8 (about 12 pages) on exactly these operator estimates. *[9 Nov: LM26 prove the asymptotically flat
  analogues in full (Lemmas 8.4–8.6, Prop. 7.1), with a bound on spatial primitives (Lemma 8.6, Prop. 7.2) in place of a
  Lipschitz bound on the first time derivative of the boundary data.]*
- **Imported:** the h-principle is quoted from Mengual 2022, with the caveat stated in the text. The paper relies on the same theorem.
- **Regularity:** `C^{6,γ}` is cruder than needed; the paper works with `H⁵`.

## 6. Small observations on the actual paper

- The introduction (§1.4) states the opening constant as `0 < c̄ < 1`. The single-layer proof (§3.3.4, §7.4) gives `0 < c̄ < 1/2`. The full range relies on the multilayer extension of §7.5, which the authors present "briefly". *[9 Nov: `c̄ < 1/2` is where the
  lower bound of §7.4 stops, not the exact hull condition (4.14); see the Correction at the top.]*
- Typos: "Se Appendix A" (§2.1), "smooth functon" (§8.1).
- The paper declares AI assistance for the relaxation result in Appendix A (Λ-affine functions suggested by "ChatGPT 6") and for presentation.

## 7. Bottom line

From a one-paragraph abstract, which already announced the central idea (inverting in terms of σ rather than the vorticity), the reconstruction reproduced the overall construction of the actual paper *[Corrected 9 Nov 2026: the 6 Nov text said "the central mechanism and the overall construction"; the central idea was announced in the abstract, so the agreement on it is not independent]*, including some details the abstract did not mention:
- one intermediate density with asymmetric `1/(1∓A)`-type speeds;
- the second-order time expansion;
- the same operator inverted at both orders;
- the optimal multilayer rates, as a formal remark. *[Corrected 9 Nov 2026: the multilayer supremum speeds, as a formal remark
  written after LM26 v1 appeared, so without priority; "optimal" is not established.]*

What it lacks relative to the paper is:
- the explicit second-order formulas and full operator estimates on `R`;
- the exact characterisation of the relaxation;
- the rigorous multilayer extension;
- the sharper `H⁵` regularity;
- *[9 Nov:]* mixing in every time slice (LM26 Def. 2.2), and the flat strip with free intermediate density (LM26 §3.2).

Since the result and its announcement belong to Luna-Velasco and Mengual, the reconstruction should be treated as a study document, not as an independent publication.
