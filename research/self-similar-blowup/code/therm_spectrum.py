#!/usr/bin/env python3
"""therm_spectrum.py -- spectra of the linearised similarity operators (Section 4).

Models (similarity variables, ODE profile, sector g = r^l Y_l(omega) h(rho), rho = r^2):

  reduced  u_tt - Lap u = delta e^u      L_R(h1,h2) = (-Lam h1 + h2, Lap h1 - Lam h2 - h2 + 2 h1)
  full     (d_t)(u_t - e^u) = Lap u      L_F(h1,h2) = (-Lam h1 + h2, Lap h1 - Lam h2 + h1)

with, on the sector, Lam h = l h + 2 rho h',  Lap h = 4 rho h'' + (4 l + 2 d) h'.

Checks performed
  (A) Chebyshev collocation in rho on [0,1] (no boundary conditions) for d=1..5,
      l=0..4: eigenvalues vs. the predicted lists
          reduced: 1-l-2n, -2-l-2n ;  full: 1-l-2n, -1-l-2n.
  (B) Exact rational computation on polynomials (monomial basis, sympy):
      algebraic vs geometric multiplicities -> semisimplicity of lambda = 1, 0
      and the Jordan block of the full model at lambda = -1 (l = 0).
  (C) Independent check of the connection coefficient: integrate the radial
      ODE numerically from rho=0 (regular Frobenius solution) to rho -> 1,
      decompose in the two Frobenius solutions at rho = 1, and compare the
      coefficient of the singular one with Gamma(c)Gamma(a+b-c)/(Gamma(a)Gamma(b));
      argument-principle zero count of the numerically computed coefficient on a
      contour in Re(lambda) > -1/2.
"""
import numpy as np
import sympy as sp
import mpmath as mp


# ---------------------------------------------------------------------------
# Chebyshev collocation
# ---------------------------------------------------------------------------
def cheb(N):
    """Chebyshev differentiation matrix on x_j = cos(pi j/N), j=0..N (Trefethen)."""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def sector_matrix(d, l, model, N):
    D, x = cheb(N)
    rho = (x + 1) / 2            # map [-1,1] -> [0,1]
    D1 = 2 * D
    D2 = D1 @ D1
    I = np.eye(N + 1)
    R = np.diag(rho)
    Lam = l * I + 2 * R @ D1
    Lap = 4 * R @ D2 + (4 * l + 2 * d) * D1
    if model == 'reduced':
        M = np.block([[-Lam, I], [Lap + 2 * I, -Lam - I]])
    else:
        M = np.block([[-Lam, I], [Lap + I, -Lam]])
    return M


def predicted(model, l, lmin=-6.0):
    out = []
    for n in range(0, 20):
        a = 1 - l - 2 * n
        b = (-2 - l - 2 * n) if model == 'reduced' else (-1 - l - 2 * n)
        out += [v for v in (a, b) if v > lmin]
    return sorted(out, reverse=True)


def part_A(cut=-1.6):
    """Collocation is exact on polynomials (L maps P_N x P_N into itself), but the matrix is
    highly non-normal, so in floating point only eigenvalues with Re > cut are resolved;
    below that, rounding produces spurious complex pairs (pseudospectral pollution)."""
    print(f"(A) Chebyshev collocation in rho, eigenvalues with Re > {cut}")
    worst = 0.0
    for model in ['reduced', 'full']:
        for d in [1, 2, 3, 4, 5]:
            for l in [0, 1, 2, 3, 4]:
                if d == 1 and l > 1:
                    continue
                lamN = np.linalg.eigvals(sector_matrix(d, l, model, 32))
                e = np.sort(lamN[lamN.real > cut].real)[::-1]
                im = np.max(np.abs(lamN[lamN.real > cut].imag)) if len(e) else 0.0
                pred = predicted(model, l, cut)
                if len(pred) == 0 and len(e) == 0:
                    err = 0.0
                elif len(e) == len(pred):
                    err = float(np.max(np.abs(np.sort(lamN[lamN.real > cut]) - np.sort(pred))))
                else:
                    err = np.inf
                worst = max(worst, err)
                if d in (1, 3) and l in (0, 1, 2):
                    print(f"  {model:7s} d={d} l={l}: N=32 -> {np.round(e, 5)}  (max|Im| {im:.1e})"
                          f"   predicted {pred}   max err {err:.1e}")
    print(f"  worst deviation over d=1..5, l=0..4, both models: {worst:.2e}")
    print("  (double eigenvalues of the full model at -1 (l=0), -2 (l=1) split by O(sqrt(eps)):"
          " signature of a Jordan block)")


# ---------------------------------------------------------------------------
# (B) exact multiplicities on polynomials
# ---------------------------------------------------------------------------
def poly_matrix(d, l, model, N):
    """Exact matrix of L on P_N x P_N in the monomial basis rho^m (sympy Rationals)."""
    d = sp.Rational(d)
    n = N + 1
    M = sp.zeros(2 * n, 2 * n)
    for m in range(n):
        # images of (rho^m, 0) and (0, rho^m)
        lam_m = l + 2 * m                        # Lam rho^m = (l+2m) rho^m
        lap_c = 4 * m * (m - 1) + (4 * l + 2 * d) * m  # Lap rho^m = lap_c rho^(m-1)
        # column for (rho^m, 0)
        M[m, m] += -lam_m
        if m >= 1:
            M[n + m - 1, m] += lap_c
        M[n + m, m] += 2 if model == 'reduced' else 1
        # column for (0, rho^m)
        M[m, n + m] += 1
        M[n + m, n + m] += -lam_m - (1 if model == 'reduced' else 0)
    return M


def part_B():
    print("\n(B) exact multiplicities on P_N x P_N (sympy, N = 6)")
    for model in ['reduced', 'full']:
        for d in [1, 2, 3, 4]:
            for l in [0, 1, 2]:
                if d == 1 and l > 1:
                    continue
                M = poly_matrix(d, l, model, 6)
                info = []
                for lam in [1, 0, -1, -2, -3]:
                    A = M - lam * sp.eye(M.shape[0])
                    g = M.shape[0] - A.rank()
                    if g == 0:
                        continue
                    alg = M.shape[0] - (A**(2 * M.shape[0] // 2)).rank() if False else None
                    # algebraic multiplicity = nullity of A^k for k large (k = 4 suffices here)
                    alg = M.shape[0] - (A**4).rank()
                    info.append(f"lam={lam}: geo {g}, alg {alg}")
                if d in (1, 3) or model == 'full':
                    print(f"  {model:7s} d={d} l={l}: " + "; ".join(info))


# ---------------------------------------------------------------------------
# (C) connection coefficient, computed by direct numerical integration
# ---------------------------------------------------------------------------
mp.mp.dps = 20


def regular_solution_at(alpha, beta, c, rho1):
    """Integrate rho(1-rho)h'' + (c-(a+b+1)rho)h' - ab h = 0 from rho0 small (series start)
    to rho1 with mpmath's Taylor ODE solver; return h(rho1), h'(rho1)."""
    rho0 = mp.mpf('0.05')
    # Frobenius series at 0 (truncated 2F1 sum; converges fast at rho0 = 0.05)
    terms = 40
    h0, dh0, coef = mp.mpf(0), mp.mpf(0), mp.mpf(1)
    for k in range(terms):
        h0 += coef * rho0**k
        if k > 0:
            dh0 += k * coef * rho0**(k - 1)
        coef = coef * (alpha + k) * (beta + k) / ((c + k) * (k + 1))
    f = lambda r, Y: [Y[1], -((c - (alpha + beta + 1) * r) * Y[1] - alpha * beta * Y[0]) / (r * (1 - r))]
    sol = mp.odefun(f, rho0, [h0, dh0])
    Y = sol(rho1)
    return Y[0], Y[1]


def connection_numeric(alpha, beta, c, eps=mp.mpf('0.04')):
    rho1 = 1 - eps
    h, dh = regular_solution_at(alpha, beta, c, rho1)
    z = 1 - rho1
    e = c - alpha - beta
    # Frobenius solutions at rho = 1 (functions of z = 1 - rho), derivative w.r.t. rho = -d/dz
    ua = lambda zz: mp.hyp2f1(alpha, beta, alpha + beta - c + 1, zz)
    us = lambda zz: zz**e * mp.hyp2f1(c - alpha, c - beta, e + 1, zz)
    Ua, Us = ua(z), us(z)
    dUa, dUs = -mp.diff(ua, z), -mp.diff(us, z)
    Amat = mp.matrix([[Ua, Us], [dUa, dUs]])
    coeffs = mp.lu_solve(Amat, mp.matrix([h, dh]))
    return coeffs[1]


def B_exact(alpha, beta, c):
    return mp.gamma(c) * mp.gamma(alpha + beta - c) / (mp.gamma(alpha) * mp.gamma(beta))


def ab(model, lam, l):
    if model == 'reduced':
        return (lam + l - 1) / 2, (lam + l + 2) / 2
    return (lam + l - 1) / 2, (lam + l + 1) / 2


def part_C():
    print("\n(C) connection coefficient: direct ODE integration vs closed form")
    rng = np.random.default_rng(1)
    for model in ['reduced', 'full']:
        for (d, l) in [(2, 0), (3, 1), (5, 2)]:
            c = mp.mpf(l) + mp.mpf(d) / 2
            errs = []
            for _ in range(3):
                lam = mp.mpc(rng.uniform(-0.4, 1.6), rng.uniform(-2, 2))
                a_, b_ = ab(model, lam, l)
                Bn = connection_numeric(a_, b_, c)
                Be = B_exact(a_, b_, c)
                errs.append(float(abs(Bn - Be) / abs(Be)))
            print(f"  {model:7s} d={d} l={l}: rel. errors {['%.1e' % e for e in errs]}")
    # argument principle on the contour Re lam in [-0.45, 2.45], Im in [-3, 3]
    print("  zero count of the numerically computed B(lambda) in the box [-0.45,2.45]x[-3,3]:")
    # dimensions chosen so that c-a-b is not an integer at integer lambda (no log case on
    # the contour/at the zeros): reduced d = 2, full d = 3
    for model, d in [('reduced', 2), ('full', 3)]:
        for l in [0, 1, 2]:
            c = mp.mpf(l) + mp.mpf(d) / 2
            pts = []
            box = [(-0.45, -3), (2.45, -3), (2.45, 3), (-0.45, 3), (-0.45, -3)]
            M = 24
            for (x0, y0), (x1, y1) in zip(box[:-1], box[1:]):
                for j in range(M):
                    s_ = j / M
                    pts.append(mp.mpc(x0 + s_ * (x1 - x0), y0 + s_ * (y1 - y0)))
            pts.append(pts[0])
            vals = []
            for lam in pts:
                a_, b_ = ab(model, lam, l)
                # divide out the known factor Gamma(a+b-c) (poles of the normalisation of the
                # singular Frobenius solution); zeros of D are the eigenvalues
                vals.append(complex(connection_numeric(a_, b_, c) / mp.gamma(a_ + b_ - c)))
            vals = np.array(vals)
            winding = np.sum(np.angle(vals[1:] / vals[:-1])) / (2 * np.pi)
            print(f"    {model:7s} d={d} l={l}: winding number {winding:+.3f}  "
                  f"(predicted zeros: {[z for z in predicted(model, l, -0.45) if z < 2.45]})")


if __name__ == "__main__":
    part_A()
    part_B()
    part_C()
