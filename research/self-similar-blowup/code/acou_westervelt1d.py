#!/usr/bin/env python3
"""1D numerics for the inviscid Westervelt equation near its degeneracy.

Pressure form (1-2kp)p_tt - c^2 p_xx = 2k p_t^2 is equivalent (q = 1-2kp, c = 1) to
        (q^2)_tt = 2 q_xx,
solved here as the first-order system  G_t = V,  V_t = 2 (sqrt G)_xx,  G = q^2,
periodic on [-pi, pi], second-order finite differences, RK4, CFL dt = 0.3 dx sqrt(min q).

Data: q(0) = 1 (p = 0), q_t(0) = -B (1 + eps cos x)  (a coherent compression ramp with a
mild spatial modulation).  Formal prediction of Section 7 (linearisation at the
homogeneous degenerate solution q^2 = A(T-t), similarity variable eta = x/(T-t)^{3/4}):
near the first degeneracy point (T, 0)
        q^2(t,x) = A rho + H x^2 + C rho^{3/2} + O(rho^2, rho x^2, x^4),   rho = T - t,
with the curvature mode (lambda = -1/2) forcing   C = (8/3) d H A^{-1/2}   (d = 1, c = 1).
In particular q^2 is C^1 but not C^2 in t at the degeneracy point, p -> 1/(2k) stays
bounded while p_t ~ rho^{-1/2}.  We fit A, C (time series at x = 0) and H (spatial
profile) and compare C with (8/3) H / sqrt(A).

Run: python3 acou_westervelt1d.py
"""
import numpy as np
from scipy.optimize import curve_fit


def rhs(G, V, dx):
    q = np.sqrt(np.maximum(G, 0.0))
    qxx = (np.roll(q, -1) - 2*q + np.roll(q, 1))/dx**2
    return V, 2.0*qxx


def run(N=4096, B=1.0, eps=0.2, Gstop=2e-5):
    x = np.linspace(-np.pi, np.pi, N, endpoint=False)
    dx = x[1] - x[0]
    G = np.ones(N)
    V = -2.0*B*(1 + eps*np.cos(x))          # (q^2)_t = 2 q q_t
    t = 0.0
    i0 = np.argmin(np.abs(x))
    ts, G0, prof_t, prof = [], [], [], []
    next_snap = 0
    while G.min() > Gstop:
        dt = 0.3*dx*np.sqrt(np.sqrt(G.min()))
        k1 = rhs(G, V, dx)
        k2 = rhs(G + dt/2*k1[0], V + dt/2*k1[1], dx)
        k3 = rhs(G + dt/2*k2[0], V + dt/2*k2[1], dx)
        k4 = rhs(G + dt*k3[0], V + dt*k3[1], dx)
        G = G + dt/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
        V = V + dt/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        t += dt
        ts.append(t)
        G0.append(G[i0])
        if G[i0] < 10.0**(-1 - 0.25*next_snap) and next_snap < 12:
            prof_t.append(t)
            prof.append(G.copy())
            next_snap += 1
    return x, np.array(ts), np.array(G0), np.array(prof_t), prof


def main():
    for N in (2048, 4096, 8192):
        x, ts, G0, pt, prof = run(N=N)
        # fit q^2(t,0) = A rho + C rho^{3/2} + D rho^2 + E rho^{5/2} on 2e-4 < q^2 < 0.05
        # (exponents 1 + {0, 1/2, 1, 3/2}: lambda = -1/2 mode, its square, lambda = -3/2)
        sel = (G0 > 2e-4) & (G0 < 0.05)
        f = lambda t, T, A, C, D, E: (A*np.abs(T - t) + C*np.abs(T - t)**1.5
                                      + D*(T - t)**2 + E*np.abs(T - t)**2.5)
        p0 = (ts[-1] + G0[-1]/abs((G0[-1] - G0[-2])/(ts[-1] - ts[-2])), 1.0, 0.0, 0.0, 0.0)
        (T, A, C, D, E), cov = curve_fit(f, ts[sel], G0[sel], p0=p0, maxfev=20000)
        # also a fit without the rho^{3/2} term, for comparison of residuals
        f2 = lambda t, T, A, D, E: A*(T - t) + D*(T - t)**2 + E*(T - t)**3
        (T2, A2, D2, E2), _ = curve_fit(f2, ts[sel], G0[sel], p0=(T, A, 0, 0), maxfev=20000)
        r1 = np.max(np.abs(f(ts[sel], T, A, C, D, E) - G0[sel]))
        r2 = np.max(np.abs(f2(ts[sel], T2, A2, D2, E2) - G0[sel]))
        # H from spatial profiles: G(t,x) - G(t,0) ~ H x^2 for |x| <= 0.05
        Hs = []
        for tt, Gp in zip(pt, prof):
            m = np.abs(x) <= 0.05
            c2 = np.polyfit(x[m], Gp[m] - Gp[np.argmin(np.abs(x))], 4)
            Hs.append(c2[2])
        H = Hs[-1]
        print(f"N={N}: T={T:.6f} A={A:.5f} C={C:.4f} D={D:.3f} | H(t->T)={H:.4f} "
              f"(H at earlier snapshots {np.round(Hs[-4:], 4)}) | predicted C=(8/3)H/sqrt(A)="
              f"{8/3*H/np.sqrt(A):.4f} | fit residual with half-integer powers: {r1:.1e}, integer powers only: {r2:.1e}")
        # local-ODE comparison
        print(f"        local ODE degeneracy time 1/(2B(1+eps)) = {1/(2*1.2):.6f}")


if __name__ == "__main__":
    main()
