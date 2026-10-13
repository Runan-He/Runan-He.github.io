#!/usr/bin/env python3
"""
num_mems_threshold.py -- static and (finite-horizon) dynamic pull-in values for
w_tt + gam w_t - Delta w = -lam w^{-2}, w = 1 on the boundary, w(0) = 1,
w_t(0) = 0, on (-1,1) (d = 1) and the unit disk (d = 2, radial).
  static:  max over the branch of  Delta w = lam w^{-2}, w = 1 on the boundary
           (shooting in the scaled variable W'' + (d-1)/rho W' = W^{-2}, W(0)=1);
  dynamic: bisection in lam for "touchdown (min w < 0.05) before t = TMAX"
           with a coarse leapfrog (the undamped threshold depends on TMAX in
           principle; we report the bracket).
Runtime < 1 min.  Prints a table.
"""
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from num_common import leapfrog1d, leapfrog_rad  # noqa: E402

TMAX = 30.0


def static_pullin(d):
    # W'' + (d-1)/rho W' = W^-2, W(0) = 1;  w(r) = m W(r sqrt(lam/m^3)),
    # w(1) = 1  <=>  W(rho1) = 1/m,  lam = m^3 rho1^2.
    def f(r, y):
        return [y[1], y[0] ** -2 - ((d - 1) / r * y[1] if r > 0 else 0.0)]
    r0 = 1e-6
    y0 = [1 + r0 ** 2 / (2 * d), r0 / d]
    sol = solve_ivp(f, (r0, 60), y0, dense_output=True, rtol=1e-11, atol=1e-12)
    rho = np.linspace(r0, 60, 200001)
    W = sol.sol(rho)[0]
    m = 1 / W
    lam = m ** 3 * rho ** 2
    k = np.argmax(lam)
    return lam[k], m[k]


def touches(lam, gam, d):
    F = lambda w: -lam / w ** 2
    hit = [False]

    def stop(W):
        if W.min() < 0.05:
            hit[0] = True
            return True
        return False

    tt = [0.0]

    def rec(t, W):
        tt[0] = t

    if d == 1:
        M = 800
        stopper = lambda W: stop(W) or tt[0] > TMAX
        leapfrog1d(F, gam, M, record=rec, stop=stopper)
    else:
        M = 400
        stopper = lambda W: stop(W) or tt[0] > TMAX
        leapfrog_rad(F, gam, M, d=2, record=rec, stop=stopper)
    return hit[0]


def dyn_threshold(gam, d, lo, hi, it=7):
    assert not touches(lo, gam, d) and touches(hi, gam, d)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if touches(mid, gam, d):
            hi = mid
        else:
            lo = mid
    return lo, hi


if __name__ == "__main__":
    for d in (1, 2):
        ls, ms = static_pullin(d)
        print("d=%d: static pull-in lam* = %.5f (w(0) = %.4f at the fold)" % (d, ls, ms))
        for gam in (0.0, 1.0):
            lo, hi = dyn_threshold(gam, d, 0.15 if d == 1 else 0.4, 1.0 if d == 1 else 1.5)
            print("      gam=%.0f: touchdown before t=%.0f iff lam > lam_dyn in [%.4f, %.4f]"
                  % (gam, TMAX, lo, hi))
