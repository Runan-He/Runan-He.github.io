#!/usr/bin/env python3
"""
num_exp_touchdown.py -- light-cone convergence test for u_tt - u_xx = e^u in 1D
on (-1,1), u = 0 at x = +-1, u(0) = 0, u_t(0) = 0 (symmetric data, blowup at
x = 0).  psi1 = u + 2 log(T-t), psi2 = (T-t) u_t; ODE profile (log 2, 2).
Prediction (num_spectra.py): symmetric sector eigenvalues 1 (time
translation), -1 (eigenfunction (1 - xi^2, -(1 + xi^2))), -2 (ODE energy
mode), ...  so the deviation should decay like (T-t)^1 with shape 1 - xi^2.
Output: figures/num_exp_touchdown.pdf + stdout.  Runtime ~ 1-3 min.
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

sys.dont_write_bytecode = True   # keep the project free of __pycache__
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from num_common import Model, lagr, richardson1d, leapfrog1d, sim1d, bisect_T  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
NS, RTOL = 20, 1e-12


def main():
    t0 = time.time()
    m = Model("exp")
    rec = []
    leapfrog1d(m.F, 0.0, 2000, bc=0.0, record=lambda t, U: rec.append((t, U.max())),
               stop=lambda U: U.max() > 4.0)
    t1 = np.floor(rec[-1][0] / 0.01) * 0.01
    x, U, V, info = richardson1d(m.F, 0.0, t1, bc=0.0)
    print("u_tt - u_xx = e^u, (-1,1), u=0 at the ends, zero data;  t1 = %.2f, u(t1,0) = %.4f, "
          "Richardson error %.1e" % (t1, U[len(x) // 2], info["dW_rich"]))
    uf = lambda X: lagr(x, U, X)
    vf = lambda X: lagr(x, V, X)
    ev = lambda t, y: y[0] - 60
    ev.terminal = True
    s = solve_ivp(lambda t, y: [y[1], np.exp(y[0])], (0, 10), [U[len(x) // 2], V[len(x) // 2]],
                  events=ev, rtol=1e-10, atol=1e-12)
    Test = t1 + s.t_events[0][0]
    run = lambda T, dense=False, Ns=NS: sim1d(m, uf, vf, T, t1, 0.0, Ns=Ns, rtol=RTOL, dense=dense)
    lo, hi = bisect_T(run, t1 + 0.6 * (Test - t1), t1 + 1.4 * (Test - t1))
    T = 0.5 * (lo + hi)
    lo2, hi2 = bisect_T(lambda T: run(T, Ns=NS + 8), lo - 1e-7, hi + 1e-7)
    Tm = T - t1
    print("T (bisection) = %.15f ; Ns+8: %.15f ; local ODE estimate %.6f ; cone radius %.3f, "
          "front at |x| = %.2f" % (T, 0.5 * (lo2 + hi2), Test, Tm, abs(1 - t1 % 2)))
    sols = [run(TT, dense=True) for TT in (lo, hi)]
    xi = sols[0][2]
    n = len(xi)
    i0 = n // 2
    taus = np.linspace(0, min(sols[0][1].t[-1], sols[1][1].t[-1]), 1000)
    Y = [s_[1].sol(taus) for s_ in sols]
    P1, P2 = Y[0][:n] - m.p1, Y[0][n:] - m.p2
    d0 = 0.5 * (Y[0][i0] + Y[1][i0]) - m.p1
    agree = abs(Y[0][i0] - Y[1][i0]) < 1e-3 * abs(d0)
    j = np.argmin(agree) if not agree.all() else len(taus)
    tt, P1, P2, d0 = taus[:j], P1[:, :j], P2[:, :j], d0[:j]
    dev = np.max(abs(P1), axis=0)
    win = tt > 3

    def resid(p):
        A, om, B = p
        return (A * np.exp(-om * tt[win]) + B * np.exp(-2 * tt[win]) - d0[win]) / abs(d0[win])

    p = least_squares(resid, [d0[win][0] * np.exp(tt[win][0]), 0.9, 0.0]).x
    late = win & (tt > 0.5 * tt[win].max())
    om_late = -np.polyfit(tt[late], np.log(abs(d0[late])), 1)[0]
    # shape of the deviation at the end of the window
    from numpy.polynomial import chebyshev as Ch
    prof = P1[:, -1] / P1[i0, -1]
    cf = np.polyfit(xi ** 2, prof, 4)
    print("valid window tau in [0, %.1f] (T-t down to %.1e); deviation max|psi1 - log 2| from %.2e to %.2e"
          % (tt[-1], Tm * np.exp(-tt[-1]), dev[0], dev[-1]))
    print("omega: fit A e^{-omega tau} + B e^{-2 tau}: %.5f ; late-window slope %.5f (prediction 1)"
          % (p[1], om_late))
    print("deviation shape at end: (psi1-log2)/(psi1(0)-log2) = 1 %+.5f xi^2 %+.1e xi^4 (prediction 1 - xi^2)"
          % (cf[-2], cf[-3]))
    print("psi2 component: (psi2(xi)-2)/(psi1(0)-log2) at xi=0: %.5f, at xi=1: %.5f "
          "(prediction from the eigenfunction (1-xi^2, -(1+xi^2)): -1 and -2)"
          % (P2[i0, -1] / P1[i0, -1], P2[0, -1] / P1[i0, -1]))
    make_figure(tt, Tm, dev, d0, P1, xi, i0, np.max(P1, axis=0) - np.min(P1, axis=0))
    print("total time %.0f s" % (time.time() - t0))


def make_figure(tt, Tm, dev, d0, P1, xi, i0, osc):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9})
    fig, axs = plt.subplots(1, 2, figsize=(8.5, 3.3))
    s = Tm * np.exp(-tt)
    axs[0].loglog(s, dev, label=r"$\|\psi_1-\log 2\|_\infty$")
    axs[0].loglog(s, osc, "--", label=r"$\mathrm{osc}_\xi\,\psi_1$")
    ss = np.logspace(np.log10(s.min()), np.log10(s.max()), 10)
    axs[0].loglog(ss, 0.5 * ss, "k:", lw=1, label=r"slope 1")
    axs[0].invert_xaxis()
    axs[0].set_xlabel(r"$T-t$")
    axs[0].legend(fontsize=7)
    axs[0].grid(alpha=0.3, which="both", lw=0.3)
    xx = np.linspace(-1, 1, 201)
    from numpy.polynomial import chebyshev as Ch
    for k in np.linspace(0, len(tt) - 1, 5).astype(int)[1:]:
        cf = Ch.chebfit(xi, P1[:, k] / P1[i0, k], len(xi) - 1)
        axs[1].plot(xx, Ch.chebval(xx, cf), label=r"$T-t=%.0e$" % s[k])
    axs[1].plot(xx, 1 - xx ** 2, "k:", lw=1.2, label=r"$1-\xi^2$")
    axs[1].set_xlabel(r"$\xi$")
    axs[1].set_ylabel(r"$(\psi_1-\log2)/(\psi_1(\cdot,0)-\log 2)$")
    axs[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "num_exp_touchdown.pdf"))
    plt.close(fig)
    print("wrote figures/num_exp_touchdown.pdf")


if __name__ == "__main__":
    main()
