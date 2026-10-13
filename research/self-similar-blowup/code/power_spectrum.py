#!/usr/bin/env python3
"""Chebyshev collocation check of the eigenvalue formula (Section 5).

We discretise the *first-order* linearised operator L = L0 + L' around the ODE
profile directly (no hypergeometric reduction is used), in the angular sector
psi_j = r^l Y_l(omega) f_j(z), z = r^2 in [0,1]:

  (L f)_1 = -2 z f1' - (s+l) f1 + f2
  (L f)_2 = 4 z f1'' + 2(d+2l) f1' - 2 z f2' - (s+1+l) f2 + V0 f1

where V0 = F'(kappa) and kappa > 0 solves s(s+1) kappa = F(kappa); both are
computed *numerically* from F (so the sign conventions for the MEMS force
F(w) = -w^{-q} are tested, not assumed).

Predicted:  lambda = 1 - l - 2n  and  lambda = -2s - 2 - l - 2n,  n = 0,1,2,...

Remark (why plain collocation is "too good"): L maps P_N x P_N (polynomials of
degree <= N in z) into itself, and on the top-degree coefficients it acts as the
2x2 matrix [[-s-l-2j, 1],[V0, -s-1-l-2j]], whose eigenvalues are exactly
1-l-2j and -2s-2-l-2j.  Hence in exact arithmetic polynomial collocation returns
precisely the predicted list and can never detect a non-polynomial eigenfunction;
in floating point the lower eigenvalues are extremely ill-conditioned (spurious
"continuum" for Re(lambda) < (d-1)/2 - s, where the singular branch (1-z)^mu,
mu = (d-1)/2 - s - lambda, is continuous).  Collocation is therefore used here
only as a sign/V0 sanity check (well-conditioned eigenvalues 1, 0, -1, -2s-2).

The independent test of the eigenvalue list (completeness!) is a MONODROMY /
Evans-function test: the solution regular at z=0 (Frobenius series, recursion
derived from the ODE) is continued numerically (DOP853, rtol 1e-12) along the
circle |z-1| = 1/2 around the singular point z=1.  It returns to itself iff it
has no (1-z)^mu or log(1-z) component, i.e. iff it is analytic at z=1, i.e.
iff lambda is an eigenvalue -- except for single-valued solutions with a POLE at
z=1 (possible only when mu is a negative integer); these are detected by
computing the Laurent coefficients c_{-1}, c_{-2}, c_{-3} by contour integration
and rejected.  We scan real lambda in [-4.6, 1.3] and a box of complex lambda.

Usage: python3 power_spectrum.py      (prints tables; writes nothing)
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import solve_ivp


def cheb(N):
    """Trefethen's Chebyshev differentiation matrix on x_j = cos(pi j/N)."""
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2, np.ones(N-1), 2])*(-1)**np.arange(N+1)
    X = np.tile(x, (N+1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1/c)/(dX + np.eye(N+1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def profile_data(F, dF, p):
    s = 2.0/(p-1.0)
    g = lambda k: F(k) - s*(s+1)*k
    # bracket a positive root
    ks = np.logspace(-3, 3, 4000)
    vals = np.array([g(k) for k in ks])
    idx = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    kap = brentq(g, ks[idx[0]], ks[idx[0]+1], xtol=1e-15, rtol=1e-15)
    return s, kap, dF(kap)


def spectrum(d, l, s, V0, N):
    D, x = cheb(N)
    z = (x + 1)/2          # z in [0,1]
    Dz = 2*D
    Dz2 = Dz @ Dz
    Z = np.diag(z)
    I = np.eye(N+1)
    A11 = -2*Z @ Dz - (s+l)*I
    A12 = I
    A21 = 4*Z @ Dz2 + 2*(d+2*l)*Dz + V0*I
    A22 = -2*Z @ Dz - (s+1+l)*I
    A = np.block([[A11, A12], [A21, A22]])
    ev = np.linalg.eigvals(A)
    return ev


def converged(d, l, s, V0, N1=48, N2=72, tol=1e-6, lo=-4.6):
    e1 = spectrum(d, l, s, V0, N1)
    e2 = spectrum(d, l, s, V0, N2)
    out = []
    for e in e2:
        if e.real < lo or abs(e.imag) > 1e-6:
            continue
        if np.min(np.abs(e1 - e)) < tol:
            out.append(e.real)
    out = np.unique(np.round(out, 8))
    return out[::-1]


def predicted(l, s, lo=-4.6):
    vals = []
    for n in range(10):
        for v in (1 - l - 2*n, -2*s - 2 - l - 2*n):
            if v >= lo:
                vals.append(v)
    return np.array(sorted(set(np.round(vals, 8)), reverse=True))


CASES = {
    'p=3  (cubic, d=3)':      (3.0, lambda u: u**3, lambda u: 3*u**2),
    'p=7  (septic, d=3)':     (7.0, lambda u: u**7, lambda u: 7*u**6),
    'p=-2 (MEMS q=2, d=3)':   (-2.0, lambda w: -w**(-2.0), lambda w: 2*w**(-3.0)),
    'p=-3 (MEMS q=3, d=3)':   (-3.0, lambda w: -w**(-3.0), lambda w: 3*w**(-4.0)),
}


def regular_at0(lam, d, l, s, V0, z0=0.5, K=120):
    """Frobenius (analytic) solution at z=0 of the reduced eigen-ODE
       4z(1-z) f'' + [2(d+2l) - (4A+6) z] f' - [A(A+1) - V0] f = 0,  A = lam+s+l,
    obtained from L f = lam f by eliminating f2 = (lam+s+l) f1 + 2 z f1'."""
    A = lam + s + l
    Q = A*(A + 1) - V0
    c = np.zeros(K + 1, dtype=complex)
    c[0] = 1.0
    for k in range(K):
        c[k+1] = (4*k*(k-1) + (4*A + 6)*k + Q)*c[k]/((k + 1)*(4*k + 2*(d + 2*l)))
    k = np.arange(K + 1)
    return np.sum(c*z0**k), np.sum(c[1:]*k[1:]*z0**(k[1:] - 1))


def monodromy_defect(lam, d, l, s, V0, return_poles=False):
    """Relative change of (h, h') after continuation around |z-1| = 1/2.
    With return_poles=True also returns max_j |c_{-j}|/|h(1/2)|, j = 1,2,3, where
    c_{-j} = (1/2 pi i) oint h(z) (z-1)^{j-1} dz are the Laurent coefficients of a
    single-valued solution at z = 1 (non-zero <=> pole <=> NOT an eigenfunction)."""
    A = lam + s + l
    Q = A*(A + 1) - V0

    def rhs(th, Y):
        z = 1 - 0.5*np.exp(1j*th)
        dz = -0.5j*np.exp(1j*th)
        h, dh = Y[0], Y[1]
        d2h = (Q*h - (2*(d + 2*l) - (4*A + 6)*z)*dh)/(4*z*(1 - z))
        return [dz*dh, dz*d2h] + [dz*h*(z - 1)**(j - 1) for j in (1, 2, 3)]
    h, dh = regular_at0(lam, d, l, s, V0)
    sol = solve_ivp(rhs, [0, 2*np.pi], [complex(h), complex(dh), 0j, 0j, 0j],
                    method='DOP853', rtol=1e-12, atol=1e-14)
    Y = sol.y[:, -1]
    defect = np.hypot(abs(Y[0] - h), abs(Y[1] - dh))/np.hypot(abs(h), abs(dh))
    if return_poles:
        return defect, max(abs(Y[2:]))/(2*np.pi*abs(h))
    return defect


def real_zeros(d, l, s, V0, lo=-4.6, hi=1.3, step=0.01):
    grid = np.arange(lo, hi + step/2, step)
    vals = np.array([monodromy_defect(x, d, l, s, V0) for x in grid])
    zeros, nonzeros = [], []
    for i in range(1, len(grid) - 1):
        if vals[i] <= vals[i-1] and vals[i] <= vals[i+1]:
            r = minimize_scalar(lambda x: monodromy_defect(x, d, l, s, V0),
                                bounds=(grid[i-1], grid[i+1]), method='bounded',
                                options={'xatol': 1e-12})
            if r.fun < 1e-6:
                _, pole = monodromy_defect(r.x, d, l, s, V0, return_poles=True)
                if pole > 1e-6:      # single-valued but with a pole at z=1
                    nonzeros.append((r.x, f'pole {pole:.1e}'))
                    continue
                zeros.append((r.x, r.fun))
            else:
                nonzeros.append((r.x, f'{r.fun:.1e}'))
    return zeros, nonzeros, vals.min()


def complex_min(d, l, s, V0, xs=np.arange(-4.5, 1.31, 0.2), ys=np.arange(0.2, 3.01, 0.2)):
    m = np.inf
    for x in xs:
        for y in ys:
            m = min(m, monodromy_defect(x + 1j*y, d, l, s, V0))
    return m


def main(complex_scan=True):
    d = 3
    for name, (p, F, dF) in CASES.items():
        s, kap, V0 = profile_data(F, dF, p)
        print(f"\n=== {name}: s = {s:.6f}, kappa = {kap:.10f}, V0 = F'(kappa) = {V0:.10f}, "
              f"(s+1)(s+2) = {(s+1)*(s+2):.10f}")
        for l in range(3):
            e = spectrum(d, l, s, V0, 24)
            e = np.sort(e[(e.real > -1.2) & (abs(e.imag) < 1e-8)].real)[::-1]
            zeros, nonzeros, _ = real_zeros(d, l, s, V0)
            pre = predicted(l, s)
            print(f"  l={l}: collocation (N=24, Re>-1.2)   {np.array2string(e, precision=8)}")
            print(f"        monodromy zeros in [-4.6,1.3]  "
                  f"{np.array2string(np.array([z[0] for z in zeros])[::-1], precision=8)}"
                  f"   (max residual {max([z[1] for z in zeros] + [0]):.1e})")
            print(f"        predicted                      {np.array2string(pre, precision=8)}")
            if nonzeros:
                print(f"        (rejected local minima: {[(round(a, 4), b) for a, b in nonzeros]})")
            if complex_scan and l <= 1:
                print(f"        min defect over complex box Re in [-4.5,1.3], Im in [0.2,3]: "
                      f"{complex_min(d, l, s, V0):.3e}")
        print(f"  spectral gap omega0 = min(1, 2s+2) = {min(1.0, 2*s + 2):.6f}")
    # q = 3 (s = -1/2): the eigenvalue -1 (l=0) is a double zero (Jordan block)
    p, F, dF = CASES['p=-3 (MEMS q=3, d=3)']
    s, kap, V0 = profile_data(F, dF, p)
    print("\n  q=3, l=0: defect(-1+eps) for eps = 1e-2, 1e-3, 1e-4:",
          [f"{monodromy_defect(-1 + e, 3, 0, s, V0):.2e}" for e in (1e-2, 1e-3, 1e-4)],
          "(quadratic => double zero)")
    p, F, dF = CASES['p=-2 (MEMS q=2, d=3)']
    s, kap, V0 = profile_data(F, dF, p)
    print("  q=2, l=0: defect(-1+eps) for eps = 1e-2, 1e-3, 1e-4:",
          [f"{monodromy_defect(-1 + e, 3, 0, s, V0):.2e}" for e in (1e-2, 1e-3, 1e-4)],
          "(linear => simple zero)")


if __name__ == '__main__':
    import sys
    main(complex_scan=('--no-complex' not in sys.argv))
