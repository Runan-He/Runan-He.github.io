# Self-similar blowup in wave models from industry

Research notes (AI-assisted, multi-agent draft) on wave-type models from
engineering whose solutions form self-similar singularities. Each model's
blowup is analysed for stability with the method of Donninger et al.:
similarity coordinates, semigroup theory, spectral analysis and modulation.

**Read:** [`self-similar-blowup.pdf`](self-similar-blowup.pdf), about 74 pages.

## Contents

| Section | File | Topic |
|---|---|---|
| 1 | `sections/00_intro.tex` | Introduction, results in brief, provenance |
| 2 | `sections/01_models_literature.tex` | Survey of industrial models, Donninger-method literature, novelty check |
| 3 | `sections/02_framework.tex` | The method; universal spectrum of ODE-type blowup |
| 4 | `sections/03_mems.tex` | **Hyperbolic MEMS touchdown** (main section) |
| 5 | `sections/04_thermal.tex` | Hyperbolic (Cattaneo) thermal explosion, exponential nonlinearity |
| 6 | `sections/05_power_family.tex` | Power family, softening foundations / Klein–Gordon, damping and mass |
| 7 | `sections/06_geometric.tex` | Wave maps (director fields), relativistic membranes |
| 8 | `sections/07_acoustics.tex` | Nonlinear acoustics, derivative nonlinearities, Westervelt |
| 9 | `sections/08_numerics.tex` | Independent numerical verification |
| 10 | `sections/09_conclusions.tex` | Status table, engineering consequences, open problems |

## Main findings

- **Universal spectrum.** Take any ODE-type blowup `κ(T−t)^{−s}` with `s > −1`.
  In the light cone, its linearisation has eigenvalues `{1−ℓ−2n} ∪ {−2s−2−ℓ−2n}`,
  and the spectral gap is `min{1, 2s+2}`.
- **Hyperbolic MEMS** `w_tt + γw_t − Δw = −λw^{−q}`, `q > 1`, `d ≤ 3`.
  - Touchdown is `w ≃ c(T−t)^{2/(q+1)}` with a flat profile in the light cone,
    where `c^{q+1} = λ(q+1)²/(2(q−1))`.
  - It is stable modulo time translation and Lorentz boosts, also with damping.
  - The rate is `(T−t)^ω` for `ω < min{1, 2(q−1)/(q+1)}`, i.e. `ω < 2/3` for the
    Coulomb force.
  - The spectrum is proved; the nonlinear step is a proof sketch.
  - Focused self-similar profiles exist numerically in `d = 2, 3` and are unstable.
- **Numerics agree with the analysis.** Simulations of the clamped device give
  the exponent `0.66666`, the constant `c` to `4·10⁻⁶`, and the convergence
  rate `0.6666`. Exponents fitted at moderate gaps are biased low (0.50–0.64).

See the status table in Section 10 for the remaining models.

## Status and caveats

- The nonlinear stability arguments are **proof sketches** in the style of
  Donninger–Schörkhuber and Ostermann. They are not complete proofs.
- Network access to arXiv and journal sites was blocked while this was written.
  Several literature statements rest on abstracts or search snippets; they are
  marked `% TODO verify` in the `.tex` sources (about 20 places).
- The comparison with Gimperlein–He–Lacey (Proc. R. Soc. A 2022) uses the
  abstract only. See `sections/03_mems.tex`, the subsection "Relation to
  Gimperlein–He–Lacey".

## Rebuilding

```bash
python3 code/merge_bib.py                       # merges bib/*.bib into bib/all.bib
latexmk -pdf main.tex                           # needs texlive-latex-extra, texlive-science
```

Every script in `code/` runs from this directory with `python3 code/<name>.py`.
They need numpy, scipy, sympy, mpmath and matplotlib, and they regenerate the
figures in `figures/`. The prefixes match the sections: `mems_`, `therm_`,
`power_`, `geom_`, `acou_`, `num_`. The `*_symbolic*.py` scripts check the
formulas in their section with sympy.
