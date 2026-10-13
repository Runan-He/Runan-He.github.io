#!/usr/bin/env python3
"""
num_common.py -- shared numerical tools for the num_* scripts.

  * cheb(N)                 Chebyshev-Gauss-Lobatto differentiation matrix
  * lagr(xg, f, X)          local 10-point Lagrange interpolation, uniform grid
  * leapfrog1d(...)         CFL = 1 leapfrog for w_tt + g w_t - w_xx = F(w)
                            on (-1,1) with Dirichlet data (exact for the
                            free 1D wave equation; singularities travel
                            along grid diagonals)
  * leapfrog_rad(...)       leapfrog (CFL 1/2) for radial w_tt + g w_t - Delta w = F
                            on the unit ball of R^d, Dirichlet at r = 1
  * sim1d / simrad          evolution in similarity coordinates
                            tau = -log((T-t)/(T-t1)), xi = (x-x0)/(T-t)
                            (Chebyshev collocation, no boundary condition at
                            the light cone |xi| = 1)
  * bisect_T                bisection for the blowup time T
Not meant to be run by itself.
"""
import numpy as np
from scipy.integrate import solve_ivp


def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def lagr(xg, f, X, p=10):
    """local p-point Lagrange interpolation of grid data (uniform xg) at X."""
    h = xg[1] - xg[0]
    i0 = np.clip(np.floor((X - xg[0]) / h).astype(int) - p // 2 + 1, 0, len(xg) - p)
    s = (X - xg[i0]) / h
    out = np.zeros_like(X, dtype=float)
    for k in range(p):
        Lk = np.ones_like(X, dtype=float)
        for m in range(p):
            if m != k:
                Lk *= (s - m) / (k - m)
        out += Lk * f[i0 + k]
    return out


# ---------------------------------------------------------------------------
# physical-space solvers
# ---------------------------------------------------------------------------
def leapfrog1d(F, gam, M, t1=None, w0=None, v0=None, bc=1.0, record=None, stop=None):
    """CFL=1 leapfrog on [-1,1], h = 2/M, w = bc at x = +-1.
    F(w): nonlinearity (w_tt + gam w_t - w_xx = F(w)).
    Returns x, w(t1), w_t(t1) (if t1 given), else runs until stop(W) is True,
    recording record(t, W) every step.  w0, v0: callables (default bc, 0)."""
    h = 2.0 / M
    x = -1 + h * np.arange(M + 1)
    W0 = np.full(M + 1, bc, dtype=float) if w0 is None else w0(x)
    V0 = np.zeros(M + 1) if v0 is None else v0(x)
    W1 = W0.copy()
    W1[1:-1] = (0.5 * (W0[2:] + W0[:-2]) + 0.5 * h * (0.5 * V0[2:] + V0[1:-1] + 0.5 * V0[:-2])
                + 0.5 * h * h * (-gam * V0[1:-1] + F(W0[1:-1])))
    W1[0] = W1[-1] = bc
    a = 1 + gam * h / 2
    b = 1 - gam * h / 2
    Wm, Wc = W0, W1
    n1 = None if t1 is None else int(round(t1 / h))
    if n1 is not None:
        assert abs(n1 * h - t1) < 1e-11, "t1 must be a multiple of h"
    n = 1
    while True:
        Wn = np.empty_like(Wc)
        Wn[1:-1] = (Wc[2:] + Wc[:-2] - b * Wm[1:-1] + h * h * F(Wc[1:-1])) / a
        Wn[0] = Wn[-1] = bc
        if n1 is not None and n == n1:
            return x, Wc, (Wn - Wm) / (2 * h)
        if record is not None:
            record(n * h, Wc)
        if stop is not None and stop(Wc):
            return x, Wc, (Wn - Wm) / (2 * h)
        Wm, Wc = Wc, Wn
        n += 1


def richardson1d(F, gam, t1, Ms=(4000, 8000, 16000), **kw):
    """leapfrog1d on three grids + Richardson extrapolation (O(h^2) error)."""
    R = [leapfrog1d(F, gam, M, t1, **kw) for M in Ms]
    x = R[0][0]
    W = [r[1][:: 2 ** i] for i, r in enumerate(R)]
    V = [r[2][:: 2 ** i] for i, r in enumerate(R)]
    W1, W2 = (4 * W[1] - W[0]) / 3, (4 * W[2] - W[1]) / 3
    V1, V2 = (4 * V[1] - V[0]) / 3, (4 * V[2] - V[1]) / 3
    return x, W2, V2, dict(dW_rich=np.max(abs(W2 - W1)), dW_raw=np.max(abs(W[2] - W[1])))


def leapfrog_rad(F, gam, M, d=2, t1=None, cfl=0.5, record=None, stop=None, bc=1.0):
    """radial leapfrog, r in [0,1], h = 1/M, dt = cfl*h (conservative stencil)."""
    h = 1.0 / M
    dt = cfl * h
    r = h * np.arange(M + 1)
    rp = (r[1:-1] + h / 2) ** (d - 1)
    rm = (r[1:-1] - h / 2) ** (d - 1)
    rc = r[1:-1] ** (d - 1) * h * h

    def lap(W):
        L = np.empty_like(W)
        L[1:-1] = (rp * (W[2:] - W[1:-1]) - rm * (W[1:-1] - W[:-2])) / rc
        L[0] = 2 * d * (W[1] - W[0]) / h ** 2
        L[-1] = 0
        return L

    W0 = np.full(M + 1, bc, dtype=float)
    W1 = W0 + 0.5 * dt * dt * (lap(W0) + F(W0))
    W1[-1] = bc
    a = 1 + gam * dt / 2
    b = 1 - gam * dt / 2
    Wm, Wc = W0, W1
    n1 = None if t1 is None else int(round(t1 / dt))
    if n1 is not None:
        assert abs(n1 * dt - t1) < 1e-11
    n = 1
    while True:
        Wn = (2 * Wc - b * Wm + dt * dt * (lap(Wc) + F(Wc))) / a
        Wn[-1] = bc
        if n1 is not None and n == n1:
            return r, Wc, (Wn - Wm) / (2 * dt)
        if record is not None:
            record(n * dt, Wc)
        if stop is not None and stop(Wc):
            return r, Wc, (Wn - Wm) / (2 * dt)
        Wm, Wc = Wc, Wn
        n += 1


def even_fit(r, f, R0, deg=40):
    """least-squares fit by even Chebyshev polynomials on [0,R0] (filters the
    grid-scale noise that the 2D focusing singularity leaves in w_t)."""
    m = r <= R0
    B = np.polynomial.chebyshev.chebvander(r[m] / R0, deg)[:, ::2]
    cf = np.linalg.lstsq(B, f[m], rcond=None)[0]
    resid = np.max(abs(B @ cf - f[m]))
    return (lambda rr: np.polynomial.chebyshev.chebvander(rr / R0, deg)[:, ::2] @ cf), resid


# ---------------------------------------------------------------------------
# similarity-coordinate solvers
# ---------------------------------------------------------------------------
class Model:
    """MEMS: F = -lam w^{-q} (psi1 = (T-t)^{-sig} w),  or exp: F = e^u
    (psi1 = u + 2 log(T-t)).  Lower-order term -gam w_t."""

    def __init__(self, kind, lam=1.0, q=2, gam=0.0):
        self.kind, self.lam, self.q, self.gam = kind, lam, q, gam
        if kind == "mems":
            self.sig = 2.0 / (q + 1)
            self.c = (lam / (self.sig * (1 - self.sig))) ** (1.0 / (q + 1))
            self.p1, self.p2 = self.c, -self.sig * self.c   # ODE profile
        else:
            self.p1, self.p2 = np.log(2.0), 2.0

    def F(self, w):
        return -self.lam * w ** (-self.q) if self.kind == "mems" else np.exp(w)

    def to_sim(self, Tm, w, v):
        if self.kind == "mems":
            return Tm ** (-self.sig) * w, Tm ** (1 - self.sig) * v
        return w + 2 * np.log(Tm), Tm * v

    def rhs_terms(self, a, b):
        """(zeroth-order part of d psi1/dtau, of d psi2/dtau) without transport."""
        if self.kind == "mems":
            return self.sig * a + b, -(1 - self.sig) * b - self.lam * a ** (-self.q)
        return b - 2.0, -b + np.exp(a)

    def events(self, n):
        """continuous event functions with the flag they signal:
        +1: blowup inside the cone (T too large); -1: T too small."""
        if self.kind == "mems":
            e1 = lambda t, y: np.min(y[:n]) - 0.2 * self.p1   # quenching in the cone
            e2 = lambda t, y: np.max(y[:n]) - 5.0 * self.p1   # psi1 -> infinity
        else:
            e1 = lambda t, y: np.max(y[:n]) - (self.p1 + 3.0)  # blowup in the cone
            e2 = lambda t, y: np.min(y[:n]) - (self.p1 - 3.0)  # psi1 -> -infinity
        e1.terminal = e2.terminal = True
        return [e1, e2], [+1, -1]

    def breakdown_flag(self, a):
        if self.kind == "mems":
            return +1 if np.min(a) < self.p1 else -1
        return +1 if np.max(a) > self.p1 else -1


def _run(rhs, y0, n, model, taumax, dense, rtol):
    evs, flags = model.events(n)
    sol = solve_ivp(rhs, (0, taumax), y0, method="DOP853", rtol=rtol, atol=1e-13,
                    events=evs, dense_output=dense)
    f = 0
    for k, fl in enumerate(flags):
        if sol.t_events[k].size:
            f = fl
    if f == 0 and sol.status < 0:     # integrator broke down near the singularity
        f = model.breakdown_flag(sol.y[:n, -1])
    return f, sol


def sim1d(model, wfun, vfun, T, t1, x0, Ns=24, taumax=45.0, dense=False, rtol=1e-11):
    """1D light-cone evolution; wfun, vfun give w(t1,.), w_t(t1,.)."""
    D, xi = cheb(Ns)
    D2 = D @ D
    Tm = T - t1
    X = x0 + Tm * xi
    p1, p2 = model.to_sim(Tm, wfun(X), vfun(X))
    n = Ns + 1
    g = model.gam

    def rhs(tau, y):
        a, b = y[:n], y[n:]
        za, zb = model.rhs_terms(a, b)
        da = -xi * (D @ a) + za
        db = D2 @ a - xi * (D @ b) + zb - g * Tm * np.exp(-tau) * b
        return np.concatenate([da, db])

    f, sol = _run(rhs, np.concatenate([p1, p2]), n, model, taumax, dense, rtol)
    return f, sol, xi


def simrad(model, wfun, vfun, T, t1, d=2, Ns=16, taumax=45.0, dense=False, rtol=1e-11):
    """radial light-cone evolution centred at the origin, variable z = rho^2."""
    D, x = cheb(Ns)
    z = (1 + x) / 2
    Dz = 2 * D
    Dzz = Dz @ Dz
    Tm = T - t1
    rr = Tm * np.sqrt(z)
    p1, p2 = model.to_sim(Tm, wfun(rr), vfun(rr))
    n = Ns + 1
    g = model.gam

    def rhs(tau, y):
        a, b = y[:n], y[n:]
        az = Dz @ a
        za, zb = model.rhs_terms(a, b)
        da = -2 * z * az + za
        db = 4 * z * (Dzz @ a) + 2 * d * az - 2 * z * (Dz @ b) + zb - g * Tm * np.exp(-tau) * b
        return np.concatenate([da, db])

    f, sol = _run(rhs, np.concatenate([p1, p2]), n, model, taumax, dense, rtol)
    return f, sol, z


def bisect_T(run, lo, hi, rtol=4e-16, maxit=80):
    """run(T) -> (flag, ...); flag(lo) = -1, flag(hi) = +1."""
    flo, fhi = run(lo)[0], run(hi)[0]
    if not (flo == -1 and fhi == +1):
        raise RuntimeError("bad bracket: flags %s %s at %r %r" % (flo, fhi, lo, hi))
    it = 0
    while hi - lo > rtol * abs(hi) and it < maxit:
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:
            break
        f = run(mid)[0]
        if f == +1:
            hi = mid
        elif f == -1:
            lo = mid
        else:
            break
        it += 1
    return lo, hi
