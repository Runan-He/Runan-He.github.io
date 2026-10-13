"""geom_modes.py -- spectral problem for self-similar co-rotational wave maps.

Shared numerical kernel used by geom_wavemaps_spectrum.py.

Mode equation (radial, in R^n, "s = 1" convention of NOTATION.md, y = rho^2):
    y(1-y) G'' + (n/2 - (lam+5/2) y) G' - ((lam+1)(lam+2) + V(y))/4 G = 0,
with V(y) = -A/(B+y)^p.
  original problem (wave map, w = u/r):  n = d+2, A = 8(d-1)(d-2), B = d-2, p = 2
  SUSY partner (gauge mode removed):     n = d+4, A = 4(d-1),      B = d-2, p = 1,
                                          spectral parameter lam' = lam + 1.
Eigenvalue  <=>  the solution analytic at y=0 is analytic (C^infty) at y=1.
"""
import numpy as np
import scipy.linalg as sl
import mpmath as mp


def problem(d, which="orig"):
    if which == "orig":
        return dict(n=d + 2, A=8 * (d - 1) * (d - 2), B=d - 2, p=2, shift=0)
    if which == "susy":
        return dict(n=d + 4, A=4 * (d - 1), B=d - 2, p=1, shift=1)
    raise ValueError(which)


# ---------------------------------------------------------------- collocation
def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(1))
    return D, x


def collocation_spectrum(d, N, which="orig"):
    """Eigenvalues lam of the quadratic pencil, via companion linearisation.
    No boundary conditions are imposed: y=0 and y=1 are both singular points and
    collocation on Chebyshev-Lobatto nodes selects solutions analytic there."""
    P = problem(d, which)
    n, A, B, p, sh = P["n"], P["A"], P["B"], P["p"], P["shift"]
    D, x = cheb(N)
    y = (1 - x) / 2
    Dy = -2 * D
    D2 = Dy @ Dy
    V = -A / (B + y) ** p
    I = np.eye(N + 1)
    # 4*(eq):  lam'^2 G + lam' (4y G' + 3G) = 4y(1-y)G'' + (2n-10y)G' - (2+V)G
    Aop = 4 * np.diag(y * (1 - y)) @ D2 + np.diag(2 * n - 10 * y) @ Dy - np.diag(2 + V)
    Bop = 4 * np.diag(y) @ Dy + 3 * I
    M = np.block([[0 * I, I], [Aop, -Bop]])
    ev = sl.eigvals(M) - sh
    return ev[np.argsort(-ev.real)]


# ---------------------------------------------------------------- Frobenius / Wronskian
def _polys(d, lam, which):
    P = problem(d, which)
    n, A, B, p = P["n"], mp.mpf(P["A"]), mp.mpf(P["B"]), P["p"]
    lam = lam + P["shift"]

    def mul(a, b):
        r = [mp.mpf(0)] * (len(a) + len(b) - 1)
        for i, u in enumerate(a):
            for j, v in enumerate(b):
                r[i + j] += u * v
        return r

    Bp = [mp.mpf(1)]
    for _ in range(p):
        Bp = mul(Bp, [B, mp.mpf(1)])
    P2 = mul([4 * c for c in Bp], [0, 1, -1])
    P1 = mul([4 * c for c in Bp], [mp.mpf(n) / 2, -(lam + mp.mpf(5) / 2)])
    P0 = [-(lam + 1) * (lam + 2) * c for c in Bp]
    P0[0] += A
    return P2, P1, P0


def _shift(poly, y0):
    r = [mp.mpf(0)] * len(poly)
    for i, a in enumerate(poly):
        for j in range(i + 1):
            r[j] += a * mp.binomial(i, j) * y0 ** (i - j)
    return r


def _frob(P2, P1, P0, y0, J):
    """Coefficients of the exponent-0 Frobenius solution about y0 (a_0 = 1)."""
    p2 = _shift(P2, y0) + [0] * 6
    p1 = _shift(P1, y0) + [0] * 6
    p0 = _shift(P0, y0) + [0] * 6
    a = [mp.mpf(1)]
    for j in range(J):
        rhs = 0
        for i in range(2, len(p2)):
            k = j - i + 2
            if 0 <= k < len(a):
                rhs += p2[i] * k * (k - 1) * a[k]
        for i in range(1, len(p1)):
            k = j - i + 1
            if 0 <= k < len(a):
                rhs += p1[i] * k * a[k]
        for i in range(len(p0)):
            k = j - i
            if 0 <= k < len(a):
                rhs += p0[i] * a[k]
        a.append(-rhs / ((j + 1) * (p2[1] * j + p1[0])))
    s2 = 1 - p1[0] / p2[1]  # second Frobenius exponent
    return a, s2


def _ev(a, z):
    v = mp.mpf(0)
    dv = mp.mpf(0)
    for j in range(len(a) - 1, -1, -1):
        v = v * z + a[j]
    for j in range(len(a) - 1, 0, -1):
        dv = dv * z + j * a[j]
    return v, dv


def wronskian(d, lam, which="orig", J=400, ym=mp.mpf(1) / 2):
    """Entire function of lam whose zeros are exactly the eigenvalues:
    W(lam) * rgamma(1 - sigma2), sigma2 = second exponent at y = 1.
    (The factor 1/Gamma cancels the poles of the y=1 series at resonances
    sigma2 in N; log-free resonances give automatic zeros, as they should.)"""
    try:
        P2, P1, P0 = _polys(d, lam, which)
        a0, _ = _frob(P2, P1, P0, mp.mpf(0), J)
        a1, s2 = _frob(P2, P1, P0, mp.mpf(1), J)
    except ZeroDivisionError:
        # exactly at a resonance: use continuity of the entire function
        h = mp.mpf(10) ** (-mp.mp.dps // 2)
        return (wronskian(d, lam + h, which, J, ym) + wronskian(d, lam - h, which, J, ym)) / 2
    g0, dg0 = _ev(a0, ym)
    g1, dg1 = _ev(a1, ym - 1)
    return (g0 * dg1 - dg0 * g1) * mp.rgamma(1 - s2)


def refine(d, guess, which="orig", dps=30):
    with mp.workdps(dps):
        return mp.findroot(lambda l: wronskian(d, l, which), mp.mpc(guess))


def count_zeros(d, box, which="orig", npts=1600, dps=25):
    """Argument principle for wronskian() on the rectangle box=(x0,x1,y0,y1)."""
    x0, x1, y0, y1 = box
    with mp.workdps(dps):
        per = 2 * (x1 - x0) + 2 * (y1 - y0)
        pts = []
        for k in range(npts):
            s = per * k / npts
            if s < x1 - x0:
                z = complex(x0 + s, y0)
            elif s < (x1 - x0) + (y1 - y0):
                z = complex(x1, y0 + s - (x1 - x0))
            elif s < 2 * (x1 - x0) + (y1 - y0):
                z = complex(x1 - (s - (x1 - x0) - (y1 - y0)), y1)
            else:
                z = complex(x0, y1 - (s - 2 * (x1 - x0) - (y1 - y0)))
            pts.append(z)
        pts.append(pts[0])
        vals = [complex(wronskian(d, mp.mpc(z), which)) for z in pts]
        dargs = np.angle(np.array(vals[1:]) / np.array(vals[:-1]))
        return float(np.sum(dargs) / (2 * np.pi)), float(np.max(np.abs(dargs)))
