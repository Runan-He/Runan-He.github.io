# Numerics for the note "Global existence of smooth self-similar quenching profiles"

* `profile_ode.py` – the profile ODE `(1-η²)v'' + ((n-1)/η - 2η/3)v' + 2v/9 - v⁻² = 0`,
  forward integration from η = 0, the Taylor series of the analytic branch at the light cone η = 1,
  and backward/forward integration from the cone.
* `inner_phase.py n` – computes the inner constants `A1, A2` (oscillation of `V0` around `a r^{2/3}`)
  and `B1, B2` (phase of the solution of the linearised equation that is regular at the cone);
  writes `phases_n<n>.npy`.
* `shooting.py n [eta_m]` – two-point shooting: locates the smooth global profiles `v_k`, prints
  `c_k = v_k(0)`, `v_k(1)`, the far-field constant `b_k`, the number of intersections with the singular
  profile `a η^{2/3}` on (0,1), and compares with the formal prediction
  `c_k = exp(-(2/(3b)) (B2 - A2 + π k))`.
* `results.csv` – the table of Section 6 of the note (n = 2,…,7).

Requires `numpy` and `scipy`. Run `inner_phase.py` before `shooting.py` for each dimension.
