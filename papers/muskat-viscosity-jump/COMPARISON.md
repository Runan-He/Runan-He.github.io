# Reconstruction vs. the actual paper

**Reconstruction:** `main.pdf` in this folder (20 pp.). It was written from the seminar abstract only and committed
(`e8efc0f`) before the actual paper was opened.

**Actual paper:** P. Luna-Velasco and F. Mengual, *Mixing zone dynamics for the Saffman–Taylor instability with viscosity jump*,
arXiv:2610.34504v1 (25 Oct 2026), 45 pp.

## 1. Verdict

The reconstruction recovers the paper's architecture almost one-to-one:

- the same reformulation (pressure / Darcy field instead of the vorticity);
- the same new Muskat formulation;
- the same type of subsolution (one intermediate density, opening along the initial normal);
- the same second-order expansion in time with the same two compatibility conditions;
- the same Fredholm/energy inversion of `I − A·P`;
- the same relaxed set (shown numerically to coincide).

With the intermediate density chosen as `ρ_m = A`, the reconstruction gives exactly the paper's single-layer mixing zone and its admissible opening range.

The paper goes further in four places:
- **Exact relaxation:** it proves the relaxation is exact.
- **Multilayers:** it pushes the opening rate to the optimal range.
- **Setting:** it works in the asymptotically flat, `H^5` setting.
- **Proofs:** it gives fully explicit second-order formulas and about 12 pages of operator estimates.

On those steps the reconstruction is shorter and partly at sketch level.

## 2. Side-by-side

| Item | Reconstruction | Luna-Velasco–Mengual |
|---|---|---|
| Title | Mixing solutions for unstable Muskat data with viscosity jump (talk title) | *Mixing zone dynamics for the Saffman–Taylor instability with viscosity jump* |
| Normalisation | `ρ=±1`, `μ=1+Aρ`, `μu=−∇p−ρe₂` | `θ=±1`, `u+Aθu=−∇p−θe₂` (identical) |
| Setting | Horizontally periodic (`T×R`) | Asymptotically flat on `R` (periodic only in their Figure 1) |
| Regularity | `z₀ ∈ C^{6,γ}` (Hölder) | `z₀−(α,0) ∈ H⁵`, chord–arc |
| Saffman–Taylor function | `σ = [∂ₙp]`, per arclength | `σ = −⟦∇p⟧·∂_αz^⊥`, per unit `α` (same object, `σ_theirs = |∂_αz| σ_mine`) |
| Muskat in σ-form | `(I−2A𝒦)σ = −2n₂`, `U = −𝒦σ` | `(I−AP)σ = −2∂_αz₁`, `(∂ₜz+Q)·∂_αz^⊥ = 0`; "(Mσ) seems to be new" |
| Dual (vorticity) form | `(I+2A𝒦)γ = −2τ₂` (Remark 2.9) | `(I+AP)ϖ = −⟦θ⟧∂_αz₂`: "both problems reduce to inverting operators of the form I ± AP" |
| Why σ helps | `v = μu` has explicit curl/div; vorticity of `u` involves `curl m` | "pressure remains harmonic in each region, whereas the mixing zone generally carries interior vorticity", `(1+A²/2) curl u = AΔg` |
| Relaxed set | `{|ρ|<1, |D−a|<|a|}`; shown ⊂ second-order laminates | `O`, given by `g(z) := (1−Aθ)|m−θu|² + (1−θ²)(m−θu)·(Au+e₂) < 0`; proves `K^{lc,Λ} = K^Λ = R = O` (Appendix A, with AI-suggested Λ-affine functions) |
| h-principle | Imported from Mengual (APDE 2022) | Same source, their Theorem 2.3 |
| Subsolution | One strip, density `ρ_m ∈ (−1,1)` arbitrary; speeds `(1±ρ_m)κ` | One strip, `θ_mix = A`; speeds `c/(1∓A)` |
| Opening rate | `0 < κ < −σ₀/(2(1−Aρ_m))` | `c = −c̄ σ°/|∂_αz°|`, `c̄ < 1/2` (single layer) |
| Strip momentum | Stream function `Ψ`, `D = ∇^⊥Ψ`, linear in λ between boundary values from Rankine–Hugoniot | Potential `g`, `m = θu − (1−θ²)/(1−Aθ)(∇^⊥g + ½(Au+i))`, linear in λ between boundary values from Rankine–Hugoniot |
| Central curve | `Z = z₀ + (tU₀ + ½t²a)n₀ + h(t)e₂` | `z = z° + t z^{(1)} + ½t² z^{(2)}` |
| Order `t⁰` | Flux mismatch `M(0)=0` ⇔ Muskat velocity | `∂_αG₊ = ∂_αG₋` at `t=0` ⇔ `z^{(1)} = −Q₀[σ^{(0)}]` |
| Order `t¹` | `(I−2A𝒦₀)a = R₁` (remainder only shown to exist and be `C^{3,γ}`) | Explicit `(I−AP)σ^{(1)} = …`, `z^{(2)}` from `(5.16)` with explicit remainders `r_σ, r_G, R_a` |
| Key smallness | `Ψ⁺−Ψ⁻ = ∫M = O(t²)` | `∂_λG = ½(G₊−G₋) = O(t²)` |
| Mass / primitive | Mass identity `∫_T M = −2Q/(1−ρ_m²)` plus a vertical shift `h(t) = O(t³)` | Primitive estimates on `R` (Lemma 8.6), no shift needed |
| Inversion of `I − A·P` | Riesz–Schauder + harmonic-energy argument (Hölder) | Riesz–Schauder + harmonic-energy argument on the adjoint (`L²`, then `H^j`) |
| Multilayer / optimal rate | Formal remark (§6.2): speeds up to `−σ/(1−A)`, `−σ/(1+A)` | Proved (§7.5): `2N` layers reach any `c̄ < 1`, i.e. exactly those speeds |
| Flat case | Continuous profile with flux `θ(1−ρ²)/(1−Aρ)`, fronts `2θ/(1∓A)` | One strip: `h_± = 1/(1∓A)`; continuous profile flux `F_A = (1−θ²)/(1−Aθ)` cited from Mengual 2022 |
| Extras | Numerical/symbolic checks (`checks/`) | Vortex-blob simulation (Fig. 1); remark on spontaneous stochasticity |

## 3. What matched

1. **The key idea.** The paper says it inverts "in terms of the Saffman–Taylor function rather than the vorticity" because the pressure (equivalently my `v = μu`) is harmonic in each region. That is exactly the reconstruction's rationale. Both papers also notice the duality of the two formulations (`I ± A·P`).
2. **Same equations.** The σ-equations agree once normalisations are matched (`P ↔ 2𝒦`, `σ_theirs = |∂_αz| σ_mine`). Their boundary densities `σ_a^{(0)} = (1−aA)σ^{(0)}/2 − aAc|∂_αz°|` are the negatives of the source parts of my densities `f^±` at `t=0` with `ρ_m = A`.
3. **Same relaxation.** On 200,000 random states, my condition `D·(D−2a) < 0` and their `g < 0` agree every time; the ratio of the two quantities depends only on `(θ, A)`. Inside the strip with `θ = A` the equivalence `|D−a|<|a| ⇔ 2|∇g| < |Au+e₂|` is exact algebra.
4. **Same geometry and admissible range.** With `ρ_m = A` and `κ = c/(1−A²)`, my curves `z^± = Z ± t(1±A)κ n₀` are their `z_± = z ± t c τ^⊥/(1∓A)`. My range `(1−A²)κ < −σ₀/2` is their `c̄ < 1/2`. The positivity margins vanish at the same endpoints.
5. **Same proof skeleton.** Both papers:
   - follow the time expansion of Förster–Székelyhidi;
   - interpolate a scalar potential linearly across the strip;
   - match first order with the Muskat velocity;
   - match second order by solving the same operator equation again;
   - reduce the tangential part of the relaxed flux to a primitive of size `O(t²)`.
6. **Optimal rates.** My formal "continuous profile" remark predicted fronts at `−σ/(1−A)` and `σ/(1+A)` relative to the central curve. These are exactly the limits of their multilayer construction ("converge to the viscosity-dependent Otto profile").

## 4. Where they differ

- **Setting and regularity.**
  - The reconstruction is periodic and Hölder (`C^{6,γ}`); the paper is asymptotically flat and Sobolev (`H⁵`).
  - The periodic choice forced an extra device: the mass identity and the vertical shift `h(t)`, to keep the strip's stream function periodic.
  - On `R` the paper controls primitives directly.
- **Intermediate density.**
  - The reconstruction allows any `ρ_m`. This is a one-parameter family of asymmetries, e.g. `ρ_m → 1` gives the optimal upper speed but a vanishing lower one.
  - The paper fixes `θ_mix = A`, which maximises the flat-case product `h₊h₋`.
- **Explicitness of the second-order step.**
  - The paper computes everything explicitly: `σ^{(1)}`, `z^{(2)}`, and the remainders `R_a, r_a, r_σ, r_G`.
  - The reconstruction (Lemma 5.4) only proves that `∂ₜM(0)` is affine in the acceleration, with linear part `−(4|∂_αz₀|/(1−ρ_m²))(I−2A𝒦₀)`. It argues the remainder exists and is `C^{3,γ}` without writing it.
- **Relaxation.**
  - The reconstruction only shows its set lies in the lamination hull; this is enough for the h-principle.
  - The paper proves the full identity `K^{lc,Λ} = K^Λ = R = O`, including the bounded version `O_M`.
- **Opening rates.**
  - The reconstruction proves only the single-strip range.
  - The paper proves `c̄ < 1` via finitely many layers (§7.5, presented briefly, with details in Lemma 8.8).

## 5. Rigor of the reconstruction, compared honestly

- **At the paper's level:**
  - the σ-formulation and invertibility (Prop. 2.6–2.7; checked numerically against the vorticity formulation to `~1e−14`);
  - the laminate computation and the relaxed set;
  - the algebra of the ansatz: jump conditions, mass identity, `t=0` limits, hull at `t=0` (all checked symbolically in `checks/`).
- **Thinner:**
  - the time-regularity lemma (5.3) and the regularity of the remainder `G` (Lemma 5.4) are a structured sketch in the appendix (derivative identities plus Plemelj–Privalov), not full estimates;
  - the tubular-neighbourhood lemma (4.1) is a sketch.
  
  The paper spends §8 (about 12 pages) on exactly these operator estimates.
- **Imported:** the h-principle is quoted from Mengual 2022, with the caveat stated in the text. The paper relies on the same theorem.
- **Regularity:** `C^{6,γ}` is cruder than needed; the paper works with `H⁵`.

## 6. Small observations on the actual paper

- The introduction (§1.4) states the opening constant as `0 < c̄ < 1`. The single-layer proof (§3.3.4, §7.4) gives `0 < c̄ < 1/2`. The full range relies on the multilayer extension of §7.5, which the authors present "briefly".
- Typos: "Se Appendix A" (§2.1), "smooth functon" (§8.1).
- The paper declares AI assistance for the relaxation result in Appendix A (Λ-affine functions suggested by "ChatGPT 6") and for presentation.

## 7. Bottom line

From a one-paragraph abstract, the reconstruction reproduced the central mechanism and the overall construction of the actual paper, including some details the abstract did not mention:
- one intermediate density with asymmetric `1/(1∓A)`-type speeds;
- the second-order time expansion;
- the same operator inverted at both orders;
- the optimal multilayer rates, as a formal remark.

What it lacks relative to the paper is:
- the explicit second-order formulas and full operator estimates on `R`;
- the exact characterisation of the relaxation;
- the rigorous multilayer extension;
- the sharper `H⁵` regularity.

Since the result and its announcement belong to Luna-Velasco and Mengual, the reconstruction should be treated as a study document, not as an independent publication.
