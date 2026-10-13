#!/usr/bin/env python3
"""
num_mems_boost.py -- asymmetric data for the 1D hyperbolic MEMS equation
w_tt - w_xx = -lam w^{-2} on (-1,1), w = 1 at x = +-1,
w(0,x) = 1 + EPS x (1 - x^2), w_t(0,x) = 0.

Questions: does the light-cone profile at the touchdown point converge to a
Lorentz-boosted ODE profile  W_a(xi) = c (1-a^2)^{-1/3} (1 + a xi)^{2/3}
with a != 0, or to the flat one?

Steps
  0. sympy: the boosted family C (T - t + a x)^sigma, C^{q+1} sigma(1-sigma)
     (1-a^2) = lam, is an exact solution (any q); numerically it is a
     stationary solution of the similarity-coordinate solver.
  1. physical solution up to t1 (CFL-1 leapfrog + Richardson);
  2. for several base points x0, the blowup time T(x0) of the backward light
     cone of (T, x0) by bisection (similarity coordinates around (T, x0));
     T(x) is the blowup curve; x* = argmin T is the touchdown point;
  3. boost parameter a(tau) of the light-cone profile at each x0, from
     W(1)/W(-1) = ((1+a)/(1-a))^{2/3}; compared with T'(x0).
Output: figures/num_mems_boost.pdf + stdout.  Runtime ~ 3-6 min (4 processes).
"""
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

sys.dont_write_bytecode = True   # keep the project free of __pycache__
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from num_common import Model, lagr, richardson1d, leapfrog1d, sim1d, bisect_T  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
LAM, GAM, EPS = 0.5, 0.0, 0.25
NS, RTOL = 24, 1e-12
OFFS = (-0.12, -0.08, -0.04, 0.0, 0.04, 0.08, 0.12)


def w0(x):
    return 1 + EPS * x * (1 - x ** 2)


def symbolic_check():
    """w = C (T - t + a x)^sigma: w_tt - w_xx = (1 - a^2) W''(Z), Z = T - t + a x.
    Check  (1-a^2) W''(Z) + lam W(Z)^{-q} = 0  for q = 2, 3, 4 (Z, lam > 0, 0 < a < 1),
    and the light-cone form (T-t)^{-sigma} w = c (1-a^2)^{-sigma/2} (1 + a xi)^sigma."""
    Z, lam, b, s, xi = sp.symbols("Z lambda b s xi", positive=True)   # b = 1 - a^2, s = T - t
    out = []
    for q in (2, 3, 4):
        sig = sp.Rational(2, q + 1)
        c = (lam / (sig * (1 - sig))) ** sp.Rational(1, q + 1)
        C = c * b ** (-sig / 2)
        W = C * Z ** sig
        res = sp.simplify(b * sp.diff(W, Z, 2) + lam * W ** (-q))
        a = sp.sqrt(1 - b)
        prof = sp.simplify(s ** (-sig) * C * (s * (1 + a * xi)) ** sig
                           - c * b ** (-sig / 2) * (1 + a * xi) ** sig)
        out.append((q, res, prof))
    # chain rule part, done symbolically once
    t, x, T, a_ = sp.symbols("t x T a", real=True)
    f = sp.Function("W")
    chain = sp.simplify(sp.diff(f(T - t + a_ * x), t, 2) - sp.diff(f(T - t + a_ * x), x, 2)
                        - (1 - a_ ** 2) * sp.Subs(sp.Derivative(f(Z), Z, 2), Z, T - t + a_ * x).doit())
    return out, chain


DATA = {}


def get_data():
    if "x" not in DATA:
        m = Model("mems", lam=LAM, q=2, gam=GAM)
        rec = []
        leapfrog1d(m.F, GAM, 2000, w0=w0, record=lambda t, W: rec.append((t, W.min())),
                   stop=lambda W: W.min() < 0.3)
        t1 = np.floor(rec[-1][0] / 0.01) * 0.01
        x, W, V, info = richardson1d(m.F, GAM, t1, w0=w0)
        DATA.update(x=x, W=W, V=V, t1=t1, err=info["dW_rich"], m=m)
    return DATA


def ode_T(w, v):
    ev = lambda t, y: y[0] - 1e-6
    ev.terminal = True
    s = solve_ivp(lambda t, y: [y[1], -LAM / y[0] ** 2 - GAM * y[1]], (0, 50), [w, v],
                  events=ev, rtol=1e-10, atol=1e-12)
    return s.t_events[0][0]


def blowup_time(x0, want_traj=False):
    D = get_data()
    x, W, V, t1, m = D["x"], D["W"], D["V"], D["t1"], D["m"]
    wf = lambda X: lagr(x, W, X)
    vf = lambda X: lagr(x, V, X)
    run = lambda T, dense=False: sim1d(m, wf, vf, T, t1, x0, Ns=NS, rtol=RTOL, dense=dense)
    Test = t1 + ode_T(float(wf(np.array([x0]))[0]), float(vf(np.array([x0]))[0]))
    lo, hi = t1 + 0.5 * (Test - t1), t1 + 1.6 * (Test - t1)
    lo, hi = bisect_T(run, lo, hi)
    out = dict(x0=x0, T=0.5 * (lo + hi), t1=t1)
    if want_traj:
        sols = [run(T, dense=True) for T in (lo, hi)]
        taus = np.linspace(0, min(sols[0][1].t[-1], sols[1][1].t[-1]), 800)
        xi = sols[0][2]
        n = len(xi)
        P = [s[1].sol(taus)[:n] for s in sols]
        c = m.c

        def boost_fit(P1):
            r = (P1[0] / P1[-1]) ** 1.5          # xi = +1 is node 0, xi = -1 node n-1
            a = (r - 1) / (r + 1)
            Wa = c * (1 - a ** 2) ** (-1 / 3) * (1 + np.outer(xi, a)) ** (2 / 3)
            return a, np.max(abs(P1 - Wa), axis=0)

        a, perr = boost_fit(P[0])
        # window: the trajectories for the two ends of the T-bracket differ by
        # less than 10% of the distance to the boosted family
        agree = np.max(abs(P[0] - P[1]), axis=0) < 0.1 * perr
        j = np.argmin(agree) if not agree.all() else len(taus)
        tt = taus[:j]
        out.update(taus=tt, a=a[:j], prof_err=perr[:j],
                   flat_err=np.max(abs(P[0][:, :j] - c), axis=0), xi=xi,
                   P_end=P[0][:, j - 1], Tm=out["T"] - t1)
        w = tt > 3
        out["rate"] = -np.polyfit(tt[w], np.log(perr[:j][w]), 1)[0] if w.sum() > 10 else np.nan
    return out


def main():
    t0 = time.time()
    out, chain = symbolic_check()
    print("sympy: w_tt - w_xx - (1-a^2) W''(T-t+ax) =", chain)
    for q, res, prof in out:
        print("sympy: q=%d: residual of the boosted profile = %s ; light-cone form check = %s"
              % (q, res, prof))

    D = get_data()
    m, t1, x = D["m"], D["t1"], D["x"]
    c = m.c
    print("data: w(0,x) = 1 + %.2f x(1-x^2), lam=%.2f, gam=%.1f; t1 = %.2f; Richardson error %.1e"
          % (EPS, LAM, GAM, t1, D["err"]))
    # stationarity of the exact boosted profile in the similarity solver
    a_ex, Tm = 0.4, 0.15
    Cb = c * (1 - a_ex ** 2) ** (-1 / 3)
    f, sol, xi = sim1d(m, lambda X: Cb * (Tm + a_ex * X) ** (2 / 3),
                       lambda X: -(2 / 3) * Cb * (Tm + a_ex * X) ** (-1 / 3),
                       t1 + Tm, t1, 0.0, Ns=NS, rtol=RTOL, taumax=12, dense=True)
    Wa = Cb * (1 + a_ex * xi) ** (2 / 3)
    drift = max(np.max(abs(sol.sol(tt)[: len(xi)] - Wa)) for tt in np.linspace(0, min(12, sol.t[-1]), 50))
    print("exact boosted profile a=%.1f in the similarity solver: max drift over tau in [0,%.0f]: %.1e"
          % (a_ex, min(12, sol.t[-1]), drift))

    xc = x[np.argmin(D["W"])]
    xs = [xc + o for o in OFFS]
    with Pool(4) as pool:
        Ts = pool.map(blowup_time, xs, chunksize=1)
    X0 = np.array([o["x0"] for o in Ts])
    TT = np.array([o["T"] for o in Ts])
    p = np.polyfit(X0 - xc, TT, 4)
    dp = np.polyder(p)
    roots = np.roots(dp)
    roots = roots[np.isreal(roots)].real
    xstar = xc + roots[np.argmin(abs(roots))]
    fit_res = np.max(abs(np.polyval(p, X0 - xc) - TT))
    print("argmin w(t1,.) = %.5f ; blowup curve T(x0) (bisection):" % xc)
    for o in Ts:
        print("   x0 = %+.5f   T(x0) = %.13f   T'(x0) (poly fit) = %+.6f"
              % (o["x0"], o["T"], np.polyval(dp, o["x0"] - xc)))
    print("quartic fit residual %.1e ; x* = argmin T = %.6f ; T'' (x*) = %.4f"
          % (fit_res, xstar, np.polyval(np.polyder(dp), xstar - xc)))

    xs2 = [xstar, xstar - 0.12, xstar - 0.06, xstar + 0.06, xstar + 0.12]
    with Pool(4) as pool:
        tr = pool.map(lambda_traj, xs2, chunksize=1)
    # refit the blowup curve with all 12 values of T(x0) (degree 6)
    X0 = np.array([o["x0"] for o in Ts] + [o["x0"] for o in tr])
    TT = np.array([o["T"] for o in Ts] + [o["T"] for o in tr])
    p = np.polyfit(X0 - xc, TT, 6)
    dp = np.polyder(p)
    p4 = np.polyfit(X0 - xc, TT, 5)
    print("degree-6 fit of T(x0) on 12 points: residual %.1e ; |T' (deg 6) - T' (deg 5)| <= %.1e on the points"
          % (np.max(abs(np.polyval(p, X0 - xc) - TT)),
             np.max(abs(np.polyval(dp, X0 - xc) - np.polyval(np.polyder(p4), X0 - xc)))))
    print("\nboost parameter a of the light-cone profile (estimator from W(1)/W(-1))")
    print("   x0-x*      T(x0)            a(end of window)   T'(x0)      T-t at end   |W - W_a|   |W - c|   rate of |W-W_a|")
    for o in tr:
        Tp = np.polyval(dp, o["x0"] - xc)
        print("   %+.5f   %.13f   %+.7f         %+.7f   %.1e     %.1e     %.1e   %.3f"
              % (o["x0"] - xstar, o["T"], o["a"][-1], Tp, o["Tm"] * np.exp(-o["taus"][-1]),
                 o["prof_err"][-1], o["flat_err"][-1], o["rate"]))
        o["Tp"] = Tp
    make_figure(Ts, tr, p, dp, xc, xstar, c)
    print("total time %.0f s" % (time.time() - t0))


def lambda_traj(x0):
    return blowup_time(x0, want_traj=True)


def make_figure(Ts, tr, p, dp, xc, xstar, c):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9})
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
    X0 = np.array([o["x0"] for o in Ts] + [o["x0"] for o in tr])
    TT = np.array([o["T"] for o in Ts] + [o["T"] for o in tr])
    xx = np.linspace(X0.min() - 0.02, X0.max() + 0.02, 200)
    axs[0].plot(X0, TT, "ko", ms=4, label=r"$T(x_0)$ by bisection")
    axs[0].plot(xx, np.polyval(p, xx - xc), "k-", lw=0.8, label="polynomial fit")
    axs[0].axvline(xstar, color="r", ls="--", lw=0.8, label=r"$x^*=%.4f$" % xstar)
    axs[0].set_xlabel(r"$x_0$")
    axs[0].set_ylabel(r"blowup curve $T(x_0)$")
    axs[0].legend(fontsize=7)
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(tr)))
    for o, col in zip(sorted(tr, key=lambda o: o["x0"]), cols):
        s = o["Tm"] * np.exp(-o["taus"])
        Tp = o["Tp"]
        axs[1].semilogx(s, o["a"], color=col, label=r"$x_0-x^*=%+.2f$" % (o["x0"] - xstar))
        axs[1].axhline(Tp, color=col, ls=":", lw=1)
        axs[2].loglog(s, o["prof_err"], color=col)
    axs[1].invert_xaxis()
    axs[1].set_xlabel(r"$T(x_0)-t$")
    axs[1].set_ylabel(r"boost parameter $a(t)$ (dotted: $T'(x_0)$)")
    axs[1].legend(fontsize=6.5)
    ss = np.logspace(-8, -1, 10)
    axs[2].loglog(ss, 0.3 * ss ** (2 / 3), "k:", lw=1)
    axs[2].text(1e-6, 1e-5, r"$(T-t)^{2/3}$", fontsize=8)
    axs[2].invert_xaxis()
    axs[2].set_xlabel(r"$T(x_0)-t$")
    axs[2].set_ylabel(r"$\|W-W_{a(t)}\|_\infty$")
    for ax in axs:
        ax.grid(alpha=0.3, which="both", lw=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "num_mems_boost.pdf"))
    plt.close(fig)
    print("wrote figures/num_mems_boost.pdf")


if __name__ == "__main__":
    main()
