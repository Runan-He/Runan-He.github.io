# Novelty audit of the manuscript against arXiv:2610.34504 (LM26)

Audit date: 9 Nov 2026. Audited text: the revision of 7 Nov 2026 (commit `a1c40e3`, "Vorticity-assisted mixing
and sharp opening rates ...").

Comparison paper: P. Luna-Velasco, F. Mengual, *Mixing solutions for unstable Muskat data with viscosity jump*,
arXiv:2610.34504v1 (LM26). Line numbers `l.NNN` below refer to the plain-text extraction of the LM26 PDF that was
used for the audit (2590 lines). Equation numbers are LM26's own.

Method:
- An automated pipeline split the manuscript into 25 claims. Each claim was checked against the LM26 text by:
  - two prior-art refuters, one reading linearly and one searching by keyword;
  - a correctness reviewer;
  - a proof-breaker;
  - a judge.
- A completeness critic then looked for uncredited uses of LM26 and misdescriptions of it.
- A separate sweep covered the outside literature. It had access only to web-search snippets, because arxiv.org,
  msp.org and link.springer.com were not reachable from the session.
- The findings below that decide the verdict were re-checked by hand against the LM26 text. Two numerical checks
  were re-run independently: the LM26 range quoted in point 4 below and the git provenance.

## 1. Main finding: the paper's central contrast with LM26 is false

The revision presents a "tilting mechanism" as its main new idea. In the manuscript's account:
- LM26 use the free function in the second-order compatibility condition to make the tangential component of the
  relaxed flux D vanish;
- this is why LM26's opening range c̄ < 1/2 "ignores the vortex sheet";
- prescribing that function differently tilts D and enlarges the range to κ < κ_max.

The first two statements are wrong, and the third overstates the gain. Notation: X is the manuscript, and
ρ_m = A is LM26's intermediate density.

1. **Same state, same disc.**
   - In X's variables, LM26's mixing-zone ansatz (4.6) gives D = (Au+e₂) − 2∇^⊥g.
   - LM26's Darcy law (4.10) gives a = (e₂ − A∇p)/(1+A²).
   - Hence D·(D−2a) = (4|∇g|² − |Au+e₂|²)/(1+A²), so LM26's exact hull condition (4.14), 2|∇g| < |Au+e₂|, is X's
     disc |D−a| < |a|.
2. **LM26's layer is already tilted.**
   - LM26 choose z^(2) to cancel the O(t) part of ∂_αG₊ − ∂_αG₋ (l.712-722, Lemma 5.3). This is a condition on
     their potential g, not on X's mismatch M.
   - The result is ∂_λG = O(t²), so ∇g is tangential at leading order ((7.9), l.1390-1394), and therefore
     D·τ₀ = (Au+e₂)·τ₀ + O(t).
   - The tangential part of (4.10), together with (7.11) and (3.7b), gives
     (Au+e₂)·τ₀ = −(1−A²)γ₀/(2+A²) =: δ_LM = [2(1+A²)/(2+A²)] δ*. Here δ* = −(1−A²)γ₀/(2(1+A²)) is X's optimal
     tilt.
   - So LM26's subsolution is tilted towards the centre of the disc, and exactly onto it when A = 0.
   - On LM26's subsolution X's M is −2t∂_α(κδ_LM) + O(t²), not O(t²).
   - X's δ = 0 strip is a different subsolution from LM26's, which LM26 never use.
3. **LM26's range is lossy only in its last estimate.**
   - The opening c is kept free throughout Sections 5–7.3 (l.948-951, l.1275-1276).
   - c̄ < 1/2 enters only in §7.4 (l.1453-1457), through the lower bound |Au+e₂| ≥ (Au+e₂)·τ^⊥ = −σ₀/2. This bound
     discards the tangential (vortex-sheet) part; Prop. 7.4 computes only the normal component.
4. **What LM26's own construction already gives.**
   - With the full modulus kept and c = (1−A²)κ, LM26's unchanged single layer satisfies the strict hull at t = 0
     if and only if
     κ(2(1−A²)κ + σ₀) < (1−A²)γ₀² / (2(2+A²)²).
   - X's sharp condition at ρ_m = A is κ(2(1−A²)κ + σ₀) < (1−A²)γ₀² / (8(1+A²)).
   - The γ₀² coefficients differ by the factor 4(1+A²)/(2+A²)², which lies in (8/9, 1] and equals 1 at A = 0.
   - Independent check on X's overturned example z₀ = (α − 1.5 sin α, 0.8 sin α) (minimum over the interface):

     | A | κ_max (X) | LM26 construction, exact hull | LM26 stated range on the stable arc |
     |---|---|---|---|
     | 0 | 0.2350 | 0.2350 | empty |
     | 0.3 | 0.2297 | 0.2294 | empty |
     | −0.6 | 0.2265 | 0.2224 | empty |
     | 0.9 | 0.2959 | 0.2778 | empty |

     Pointwise, LM26's range is at least 93.8% of κ_max. Scalar example A = 0.5, σ₀ = −1, γ₀ = 1: LM26 stated
     0.6667, LM26 exact 0.7340, X 0.7347.
   - Caveat: in LM26's asymptotically flat setting, a c built from this range must also satisfy c − c_∞ ∈ H^{k−1}
     and inf c > 0. So "latent in LM26" means a few-line modification of LM26 §7.4, not a theorem LM26 state.
   - LM26 never state this range. They say their aim is "to construct mixing solutions rather than optimize the
     opening speed" (l.1482-1483).

**Consequence.** The vortex-sheet enlargement of the single-layer range, and with it mixing on Saffman–Taylor-stable
but sheared parts of the interface (pinch-free data), is latent in LM26's construction. Only LM26's final estimate
hides it. The optimal tilt δ* adds:
- a factor (2+A²)²/(4(1+A²)) ∈ [1, 9/8) on the γ₀² term;
- nothing at A = 0;
- at most about 6% in κ on the example (A = 0.9).

## 2. Other overlaps with LM26 that the revision did not credit

- **Prop. 2.10.** The relation γ = −2(τ₂ − A𝒯σ)/(1−A²) and (I+2A𝒦)γ = −2τ₂ are LM26 (3.7b), left relation,
  (1−A²)ϖ = ⟦θ⟧(AQ−i)·∂_αz, and l.455-459, in X's normalisation (ϖ = −|∂_αz|γ, ⟦θ⟧ = −2). COMPARISON.md listed
  this identity as absent from LM26, which is wrong.
- **Example 2.12 (A = 0 formulas).** σ₀ = −2n₀,₂, γ₀ = −2τ₀,₂, and "fully unstable ⟺ graph" are LM26 l.112-113 and
  l.446-447.
- **Lemma 5.2.** a^(0)·n₀ is LM26 (7.11) together with l.1432-1435. a^(0)·τ₀ is the tangential counterpart,
  obtained from (3.7b).
- **Section 5.** The O(t) mismatch is affine in the acceleration through I − AP ↔ I − 2A𝒦₀ (LM26 (5.21), (8.41)),
  and it is solvable for any right-hand side (Prop. 6.1). X's Lemmas 5.5–5.6 are LM26's second-order step with a
  nonzero target.
- **Flat single strip with arbitrary intermediate density.** LM26 §3.2 (l.573-601) has the hull |γ| < 1/2, the mass
  balance, and the choice θ_mix = A because it maximises h₊h₋. X wrote "[LM26] take ρ_m = A".
- **Relaxed set.** U_A is the set O of Mengual (Men22), written (O) in LM26. The function g in Sec. 8.3 is LM26
  (A.1).
- **Periodic operator theory.** Smoothing of 𝒦 and invertibility of I − ξ·2𝒦 for |ξ| ≤ 1 are CCG11 (LM26
  l.465-472).
- **Flat profile.** The continuous flat profile (Example 3.6) and the untilted multilayer limit (Sec. 8.1) are
  Men22 and LM26 §3.2 and §7.5. LM26 also identify the limit profile (l.1559-1564).

## 3. Provenance statement was false

The title-page note said Sections 2–4 were written from the seminar abstract before LM26 was read. git shows the
following were added in `a1c40e3` (7 Nov), after the comparison commit `164466a` (6 Nov):
- Sec. 2.5 (𝒯, Prop. 2.10, γ₀ and the pinch set);
- the constant C_Ψ and the tilt hypothesis in Lemma 4.5;
- Remark 4.6.

## 4. Outside literature (snippets only; full texts not reachable)

- **FS18 (Förster–Székelyhidi, CMP 2018).** The abstract states mixing solutions for all non-horizontal data in the
  stable regime when A = 0.
- **CFM22 (Castro–Faraco–Mengual, Ann. PDE 2022).** It treats turned interfaces at A = 0 and states "this is why we
  have opened the mixing zone slightly inside the stable regime". The confirming search did not return its exact
  admissibility inequality for the growth rate. Search snippets report |2c(α) + ∂_αz₁(α)| < 1, which in X's
  orientation is κ < (1+n₀,₂)/2 = κ_max at A = 0, and they say mixing is impossible around "fully stable points".
  These points are X's pinch set at A = 0.
- **Men22.** The abstract already says the relaxation's pinch singularity "prevents the two fluids from mixing
  wherever there is neither Rayleigh–Taylor nor vorticity at the interface".
- **No prior art found for A ≠ 0** beyond LM26. Searched for: a γ₀-dependent range, a necessity theorem, and
  curved-interface mixing on partially Saffman–Taylor-stable data.

## 5. Per-claim verdicts

| # | Claim in the revision | Verdict | Reason |
|---|---|---|---|
| C1 | Tilt mechanism (free function in the 2nd-order condition) | partly in LM26 | z^(2) and its I−AP structure are LM26's. LM26's layer is already tilted (δ_LM). New: the target treated as a design parameter, and δ* for A≠0. |
| C2 | Thm 1.1, range 0<κ<κ_max | new, minor | Not stated in LM26, but latent in their construction (§1 above). It also does not dominate LM26's multilayer range c̄<1 for fully unstable data. |
| C3 | Thm 1.2, necessity within single strips | new, minor | No necessity result for curved data in LM26. The flat case is in LM26 §3.2, and the ingredients are LM26's. |
| C4 | Cor. 1.3, pinch-free data suffice | new as a statement, minor | The existence half follows from LM26's layer with the exact hull; the "only if" half is Thm 1.2. |
| C5 | Cor. 1.4, A = 0, whole overturned interfaces | partly in LM26 (+ FS18/CFM22) | At A=0, LM26 (4.14) and (7.9) with \|e₂\|=1 give κ<(1+n₀,₂)/2 directly. |
| C6 | Cor. 1.5, injection | new, remark-level | A Galilean change of frame plus a time rescaling; it applies to LM26 Thm 1.1 just as well. G+AV is Saffman–Taylor's criterion. |
| C7 | Arbitrary ρ_m | partly in LM26 | The flat case with free θ_mix is LM26 §3.2. The curved case is routine. |
| C8 | Prop. 2.10, γ in terms of σ | already in LM26 | LM26 (3.7b) and l.455-459. |
| C9 | Relaxed disc \|D−a\|<\|a\| = O | credited, but incomplete | Men22 is not credited. |
| C10 | Periodic mass identity and vertical shift | setting difference | Correct, and specific to the periodic setting. |
| C11 | Periodic Hölder setting | setting difference | CCG11 is not credited, and the abstract dropped "horizontally periodic". |
| C12 | σ-formulation | credited | Local citations are missing. |
| C13 | Overturned example and Fig. 1 | new illustration | The framing is wrong: LM26's construction already opens on the stable arc. The caption is inaccurate, and the percentages are by parameter, not arclength. |
| C14 | Flat continuous profiles | credited (Men22) | Also LM26 §3.2 and §7.5. |
| C15 | Tilted continuous profiles (formal) | new, formal | Contains an error: the number of compatibility conditions does not grow with the number of layers (LM26 Lemma 8.8). |
| C16 | ∫σ₀ ds = −4π | setting difference | The counterpart of σ⁽⁰⁾+2 ∈ H^{k−1} in LM26. |
| C17 | Pinch set "quantifies" Men22's pinch | not established | The disc degenerates only where σ₀ = γ₀ = 0. |
| C18 | Lemma 5.2, disc centre and optimal tilt | partly in LM26 | The normal part is LM26 Prop. 7.4; the tangential part follows from (3.7b). New: the completed square and δ*. |
| C19 | Thm 1.2 applies "verbatim" to LM26 | not established | Steps 1 and 4 use periodic tools. |
| C20 | Single-strip ansatz; δ = 0 is LM26 | credited, but δ=0 ≠ LM26 | LM26's layer corresponds to δ = δ_LM. |
| C21 | I±2A𝒦 duality | credited | LM26 l.455-459; local citation missing. |
| C22 | Vortex sheet as a structural part of the relaxation | partly in LM26 | Latent in LM26 (4.14). |
| C23 | No linear opening at pinch points | new, minor | The bubble and localisation remarks overstate it. |
| C24 | Cannot be improved by enlarging the target set | new, minor | The wording was backwards. |
| C25 | Untilted profiles give LM26's multilayer limit | credited | The speeds are approached, not reached. |

## 6. What remains new relative to LM26

1. **Theorem 1.2.** Within single-strip subsolutions of constant intermediate density, the centre of mass must move
   with the Muskat velocity, and the opening must satisfy Φ(κ) ≤ 0, so κ = 0 on the pinch set. This is a short
   argument built from LM26's ingredients, but LM26 contains no such statement for curved interfaces.
2. **An explicit vortex-sheet-dependent single-layer range for A ≠ 0.** This includes openings on pinch-free,
   partially Saffman–Taylor-stable data. It is implicit in LM26's construction and not stated there or, as far as
   we found, anywhere else.
3. **The optimal tilt δ* for A ≠ 0 and general ρ_m.** It makes the range sharp within single strips, with a gain
   factor of at most 9/8 on the γ₀² term over LM26's own layer.
4. **Curved strips with ρ_m ≠ A, and the horizontally periodic Hölder setting.** Both are routine extensions.

**Overall assessment.** Relative to LM26 the manuscript's novelty is modest: a remark-level necessity theorem, plus
an explicit and slightly optimised form of a range already latent in LM26's single layer. The A = 0 statements are
very likely covered by FS18 and CFM22. As revised on 7 Nov, the manuscript was *not* more novel than LM26 in the way
its abstract claimed, and it contained false statements about LM26.

## 7. Rigour caveats (independent of novelty)

- Lemma 5.4 (time regularity) and the C^{3,γ} regularity of G are only sketched in Appendix A.4.
- The h-principle on 𝕋×ℝ is asserted, not proved.
- Mixing is shown in space–time only. LM26 prove mixing in every time slice.
- The pinch-freeness of the example for A ≠ 0 rests on numerics without validated error bounds.

## 8. Corrections to the manuscript

This file was committed first (`66d5a96`) as the audit record of the revision of 7 Nov (`a1c40e3`). The corrections based on
it are in the next commit on this branch. They:
- remove the false statements about LM26 (§1): the manuscript now says that LM26's single layer is tilted by δ_LM, that its
  exact hull condition (4.14) is the disc |D−a| < |a| and contains γ₀, that c̄ < 1/2 comes only from the last estimate of
  LM26 §7.4, and that the δ = 0 strip is a different subsolution;
- state LM26's exact-hull range κ(2(1−A²)κ+σ₀) < (1−A²)γ₀²/(2(2+A²)²) (eq. (1.8) of the manuscript, Remark 5.8) and plot it
  in Fig. 1(b) next to κ_max and LM26's stated range;
- add the missing credits of §2, in the title-page note and at each lemma, and correct the provenance note (§3);
- restate the novelty as in §6, with the new title *Sharp opening rates of single-strip subsolutions for the Muskat problem
  with viscosity jump*, in the abstract, Secs. 1–8, `COMPARISON.md` and `README.md`, and hedge the A = 0 statements and the
  FS18/CFM22 descriptions as in §4;
- fix the errors listed in §5 (total opening speed against Ott99/Sze12, arclength percentages, Fig. 1 caption, the disc on
  Π, Remark 3.3, the factor printed by `checks/compare_hulls.py`, "horizontally periodic" in the abstract, Sec. 8.1);
- add `checks/check_lm26_layer.py`, which verifies F1–F4 of this audit symbolically and on the example.

Correction to this record: the vortex-sheet subsection of `a1c40e3` is Sec. 2.5, not 2.6 (fixed in §3 above).
