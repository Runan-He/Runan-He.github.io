# Vorticity-assisted mixing and sharp opening rates for the Muskat problem with viscosity jump

This folder started as a paper reconstructed **only from the one-paragraph abstract** of the talk
*"Mixing solutions for unstable Muskat data with viscosity jump"* by Paula Luna(-Velasco), joint work with Francisco Mengual
(ICMAT, 12 November 2026); the actual paper is arXiv:2610.34504. The first version (commit `e8efc0f`) was written before that paper was read.

The current version (`main.pdf`, 24 pages) builds on arXiv:2610.34504, which it credits for the fully unstable result, and adds:
- a sharp, necessary and sufficient range for single-strip mixing zones that involves the vortex-sheet strength;
- mixing solutions for pinch-free data, including partially Saffman–Taylor stable interfaces;
- mixing along entire overturned interfaces when A = 0;
- injection-driven viscous fingering.

See `COMPARISON.md`.

* `COMPARISON.md`: what is new relative to arXiv:2610.34504, and the original comparison.
* `main.pdf`: compiled manuscript.
* `main.tex`, `sec-*.tex`, `app-cauchy.tex`, `biblio.tex`, `fig-turned.pdf`: LaTeX sources; build with `latexmk -pdf main.tex`.
* `checks/check_sigma.py`: the Saffman–Taylor and vorticity formulations give the same interface velocity (agreement ~1e-14).
* `checks/check_algebra.py`: symbolic checks of the jump conditions, mass identity, t=0 limits, flat example and laminate lemma.
* `checks/check_vorticity.py`: the vortex-strength identity and the pinch-free overturned example.
* `checks/check_sharp.py`: the sharp opening-range algebra (κ_max, optimal tilt).
* `checks/check_hull_equiv.py`: the relaxed set equals the set O of arXiv:2610.34504 (exact factor).
* `checks/compare_hulls.py`: the same equality on 200,000 random states.
* `checks/make_figure.py`: generates `fig-turned.pdf` (`python3 make_figure.py fig-turned.pdf preview.png`).
