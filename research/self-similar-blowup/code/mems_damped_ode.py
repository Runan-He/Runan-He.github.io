#!/usr/bin/env python3
"""Numerical check of the touchdown asymptotics of the (damped) MEMS ODE
        W'' + gamma W' = -lam W^{-q},   W(0) = 1, W'(0) = 0   (plate released from rest),
used in Section 3:  with theta = T - t, sigma = 2/(q+1), c^{q+1} = lam/(sigma(1-sigma)),
        W = c theta^sigma (1 + e theta^{2-2 sigma} + c1 theta + ...),   c1 = gamma/(3-q)  (q != 3),
and, for gamma = 0, e = E/(c^2 sigma (3 - 2 sigma)) with the conserved energy E = W'^2/2 - lam W^{1-q}/(q-1).
For q = 3 the exponents 1 and 2 - 2 sigma coincide and W = c theta^{1/2}(1 + (gamma/4) theta log theta + e theta + ...).
Also prints the undamped touchdown time T = int_0^1 dW / sqrt(2 lam (W^{1-q}-1)/(q-1)) (q = 2: pi/(2 sqrt(2 lam))).
The ODE is integrated with W as independent variable near touchdown (t(W), v(W)), which gives T to ~1e-12.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad


def touchdown(q, lam, gam, W_switch=0.5, W_end=1e-9):
    # phase 1: t-integration until W = W_switch
    def f(t, y):
        return [y[1], -gam * y[1] - lam * y[0] ** (-q)]

    def ev(t, y):
        return y[0] - W_switch
    ev.terminal = True
    s1 = solve_ivp(f, [0, 50], [1.0, 0.0], method="DOP853", rtol=1e-13, atol=1e-15, events=ev)
    t0, (W0, v0) = s1.t_events[0][0], s1.y_events[0][0]
    # phase 2: W as independent variable: dt/dW = 1/v, dv/dW = (-gam v - lam W^-q)/v ; integrate in log W
    def g(sv, y):
        W = np.exp(sv)
        t, v = y
        return [W / v, W * (-gam * v - lam * W ** (-q)) / v]
    Ws = np.exp(np.linspace(np.log(W0), np.log(W_end), 4001))
    s2 = solve_ivp(g, [np.log(W0), np.log(W_end)], [t0, v0], method="DOP853", rtol=1e-13, atol=1e-16,
                   t_eval=np.log(Ws))
    t = s2.y[0]
    # remaining time from W_end to 0:  W ~ c theta^sigma  ->  theta_end ~ (W_end/c)^{1/sigma}
    sig = 2.0 / (q + 1)
    c = (lam / (sig * (1 - sig))) ** (1.0 / (q + 1))
    T = t[-1] + (W_end / c) ** (1 / sig)
    return T, T - t, Ws, c, sig


def fit(theta, R, exps):
    A = np.array([theta ** e for e in exps]).T
    coef, *_ = np.linalg.lstsq(A, R, rcond=None)
    return coef


def main():
    lam = 1.0
    for q in (2, 4, 5):
        for gam in (0.0, 0.3, 1.0):
            T, theta, W, c, sig = touchdown(q, lam, gam)
            m = (theta < 2e-2) & (theta > 1e-7)
            theta, W = theta[m], W[m]
            R = W / (c * theta ** sig) - 1
            b = 2 - 2 * sig
            exps = sorted(set([round(j + k * b, 12) for j in range(3) for k in range(4)
                               if 0 < j + k * b <= 2.0 + 1e-9]))
            coef = fit(theta, R, exps)
            ce = dict(zip(exps, coef))
            line = f"q={q} gamma={gam}: T={T:.12f}  e(theta^{b:.3f})={ce[round(b, 12)]:+.6f}  " \
                   f"c1(theta^1)={ce[1.0]:+.6f}  predicted c1=gamma/(3-q)={gam / (3 - q):+.6f}"
            if gam == 0.0:
                E = -lam / (q - 1)
                line += f"   predicted e=E/(c^2 s(3-2s))={E / (c ** 2 * sig * (3 - 2 * sig)):+.6f}"
            print(line)
    # resonant case q = 3 (sigma = 1/2, 2 - 2 sigma = 1): W = c theta^{1/2}(1 + kappa theta log theta + e theta + ...),
    # kappa = gamma sigma/(2 sigma + 1) = gamma/4
    q = 3
    for gam in (0.3, 1.0):
        T, theta, W, c, sig = touchdown(q, lam, gam)
        m = (theta < 2e-2) & (theta > 1e-7)
        theta, W = theta[m], W[m]
        R = W / (c * theta ** sig) - 1
        L = np.log(theta)
        A = np.array([theta * L, theta, theta ** 2 * L ** 2, theta ** 2 * L, theta ** 2]).T
        coef, *_ = np.linalg.lstsq(A, R, rcond=None)
        print(f"q=3 gamma={gam}: kappa(theta log theta)={coef[0]:+.6f}   predicted gamma/4={gam / 4:+.6f}")
    # undamped touchdown time from rest
    for q in (2, 3, 4):
        Tq = quad(lambda W: 1 / np.sqrt(2 * lam * (W ** (1 - q) - 1) / (q - 1)), 0, 1, limit=200)[0]
        print(f"undamped touchdown time from rest, q={q}, lam=1: {Tq:.10f}" +
              (f"  (pi/(2 sqrt 2) = {np.pi / (2 * np.sqrt(2)):.10f})" if q == 2 else ""))

    # static pull-in values lambda* for Delta w = lam w^{-2} in the unit ball B_1 (d=1: interval (-1,1)),
    # w = 1 on the boundary: lam(s) = s^2 / v(s)^3 along the scaled shooting solution v'' + (d-1)/s v' = v^{-2}
    for d in (1, 2, 3):
        f = (lambda s, y, d=d: [y[1], y[0] ** -2 - ((d - 1) / s * y[1] if s > 0 else 0.0)])
        s0 = 1e-6
        sol = solve_ivp(f, [s0, 50], [1 + s0 ** 2 / (2 * d), s0 / d], rtol=1e-12, atol=1e-14, dense_output=True)
        ss = np.linspace(0.01, 50, 200000)
        v = sol.sol(ss)[0]
        lam_s = ss ** 2 / v ** 3
        print(f"static pull-in (q=2), unit ball in d={d}: lambda* = {lam_s.max():.4f};  "
              f"flat ODE touchdown from rest in the ball needs pi^2/8 = {np.pi ** 2 / 8:.4f}")


if __name__ == "__main__":
    main()
