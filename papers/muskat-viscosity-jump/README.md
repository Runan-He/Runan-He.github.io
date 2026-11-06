# Mixing solutions for unstable Muskat data with viscosity jump — reconstruction from an abstract

This folder contains a research-paper-style manuscript written **only from the one-paragraph abstract** of the talk
*"Mixing solutions for unstable Muskat data with viscosity jump"* by Paula Luna (Universidad de Sevilla), joint work with
Francisco Mengual (Seminario EDPs y Mecánica de Fluidos, ICMAT, 12 November 2026).

It was written without access to the authors' manuscript, as an exercise in reconstructing a paper from an abstract.
The announced result is due to Luna and Mengual; this text is **not** intended for submission as original work.

* `main.pdf` — compiled manuscript (20 pages).
* `main.tex` + `sec-*.tex`, `app-cauchy.tex`, `biblio.tex` — LaTeX sources (`latexmk -pdf main.tex`).
* `checks/check_sigma.py` — numerical check that the Saffman–Taylor formulation `(I-2A K)σ = -2n₂, U = -Kσ`
  gives the same interface velocity as the classical vorticity formulation `(I+2A K)γ = -2τ₂` (agreement ~1e-14).
* `checks/check_algebra.py` — symbolic checks (sympy) of the jump conditions, the mass identity, the t=0 limits,
  the relaxation inequality at t=0, the flat-interface example and the laminate lemma.
