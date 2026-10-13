#!/usr/bin/env python3
"""Search for smooth positive radial self-similar touchdown profiles of
    w_tt - Delta w = -lam w^{-q},     w = (T-t)^sigma psi(|x|/(T-t)),  sigma = 2/(q+1),
i.e. solutions of
    (1-rho^2) psi'' + ((d-1)/rho - 2(1-sigma) rho) psi' + sigma(1-sigma) psi - lam psi^{-q} = 0
on [0,1], smooth at rho = 0 and rho = 1.  By scaling psi -> lam^{1/(q+1)} psi we take lam = 1;
the constant (ODE) profile is psi = c = (1/(sigma(1-sigma)))^{1/(q+1)}.

Method (two-sided shooting):
 * left family:  psi(0) = p, psi'(0) = 0  (one parameter p > 0), integrated from rho0 to rho_m;
 * right family: the one-parameter family of solutions analytic at rho = 1 (Frobenius exponent 0;
   the other exponent is (d-1)/2 + sigma).  Parameter b = psi(1), except in the resonant case
   (d-1)/2 + sigma = 1 (i.e. d = 2, q = 3), where analyticity forces psi(1) = c and the parameter
   is b = psi'(1);
 * smooth profiles <-> intersections of the two curves  p -> (psi_L, psi_L')(rho_m),
   b -> (psi_R, psi_R')(rho_m)  in the phase plane.
Output: list of intersections for d = 1,2,3 and q = 2,3,4, and figure figures/mems_profile_shooting.pdf
"""
import os
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures", "mems_profile_shooting.pdf")
RHO_M = 0.5
ORDER = 5


def right_series(d, q):
    """Return (resonant, f) where f(param) -> list of Taylor coefficients b_j of psi in y = 1-rho."""
    sig = sp.Rational(2, q + 1)
    y = sp.symbols('y')
    bs = sp.symbols('b0:%d' % (ORDER + 2))
    psi = sum(bs[j] * y ** j for j in range(ORDER + 2))
    rho = 1 - y
    eq = (1 - rho ** 2) * sp.diff(psi, y, 2) - ((d - 1) / rho - 2 * (1 - sig) * rho) * sp.diff(psi, y) \
        + sig * (1 - sig) * psi - psi ** (-q)
    ser = sp.expand(sp.series(eq, y, 0, ORDER + 1).removeO())
    eqs = [ser.coeff(y, j) for j in range(ORDER + 1)]
    resonant = (sp.Rational(d - 1, 2) + sig == 1)
    c = (1 / (sig * (1 - sig))) ** sp.Rational(1, q + 1)
    if not resonant:
        sol = {}
        for j in range(ORDER + 1):
            e = eqs[j].subs(sol)
            sj = sp.solve(e, bs[j + 1])
            assert len(sj) == 1
            sol[bs[j + 1]] = sp.simplify(sj[0])
        funcs = [sp.lambdify(bs[0], bs[0])] + [sp.lambdify(bs[0], sol[bs[j + 1]]) for j in range(ORDER + 1)]
        return False, (lambda b: [f(b) for f in funcs])
    # resonant: eqs[0] forces b0 = c (positive root), b1 free
    sol = {bs[0]: c}
    assert sp.simplify(eqs[0].subs(sol)) == 0
    for j in range(1, ORDER + 1):
        e = eqs[j].subs(sol)
        sj = sp.solve(e, bs[j + 1])
        assert len(sj) == 1
        sol[bs[j + 1]] = sp.simplify(sj[0])
    funcs = [lambda b1: float(c), lambda b1: b1] + [sp.lambdify(bs[1], sol[bs[j + 1]]) for j in range(1, ORDER + 1)]
    return True, (lambda b: [f(b) for f in funcs])


def rhs_factory(d, q):
    sig = 2.0 / (q + 1)

    def rhs(r, Y):
        psi, dpsi = Y
        if psi <= 0:
            return [0.0, 0.0]
        dd = -(((d - 1) / r - 2 * (1 - sig) * r) * dpsi + sig * (1 - sig) * psi - psi ** (-q)) / (1 - r * r)
        return [dpsi, dd]
    return rhs


def hit_zero(r, Y):
    return Y[0] - 1e-6
hit_zero.terminal = True


def left_value(d, q, p, rho0=1e-4):
    sig = 2.0 / (q + 1)
    A = (p ** (-q) - sig * (1 - sig) * p) / (2 * d)       # psi = p + A rho^2 + ...
    Y0 = [p + A * rho0 ** 2, 2 * A * rho0]
    sol = solve_ivp(rhs_factory(d, q), [rho0, RHO_M], Y0, method="DOP853", rtol=1e-11, atol=1e-13,
                    events=hit_zero)
    if sol.status != 0 or sol.t[-1] < RHO_M - 1e-12:
        return None
    return sol.y[:, -1]


def right_value(d, q, coeffs_fun, b, eps=2e-3):
    co = coeffs_fun(b)
    y = eps
    psi = sum(cj * y ** j for j, cj in enumerate(co))
    dpsi_dy = sum(j * cj * y ** (j - 1) for j, cj in enumerate(co) if j > 0)
    Y0 = [psi, -dpsi_dy]
    if psi <= 0:
        return None
    sol = solve_ivp(rhs_factory(d, q), [1 - eps, RHO_M], Y0, method="DOP853", rtol=1e-11, atol=1e-13,
                    events=hit_zero)
    if sol.status != 0 or sol.t[-1] > RHO_M + 1e-12:
        return None
    return sol.y[:, -1]


def curve(fun, params):
    pts = []
    for s in params:
        v = fun(s)
        pts.append([np.nan, np.nan] if v is None else v)
    return np.array(pts)


def seg_intersections(P, Q, pP, pQ):
    """Intersections of polylines P (params pP) and Q (params pQ); returns list of (p, b). Vectorised."""
    out = []
    A0, A1 = P[:-1], P[1:]
    B0, B1 = Q[:-1], Q[1:]
    okA = ~(np.isnan(A0).any(1) | np.isnan(A1).any(1))
    okB = ~(np.isnan(B0).any(1) | np.isnan(B1).any(1))
    ia = np.where(okA)[0]
    ib = np.where(okB)[0]
    if len(ia) == 0 or len(ib) == 0:
        return out
    a0, da = A0[ia][:, None, :], (A1 - A0)[ia][:, None, :]
    b0, db = B0[ib][None, :, :], (B1 - B0)[ib][None, :, :]
    det = -da[..., 0] * db[..., 1] + da[..., 1] * db[..., 0]
    rhs = b0 - a0
    with np.errstate(divide='ignore', invalid='ignore'):
        s_ = (-rhs[..., 0] * db[..., 1] + rhs[..., 1] * db[..., 0]) / det
        t_ = (da[..., 0] * rhs[..., 1] - da[..., 1] * rhs[..., 0]) / det
    hit = (s_ >= 0) & (s_ <= 1) & (t_ >= 0) & (t_ <= 1) & np.isfinite(s_) & np.isfinite(t_)
    for i, j in zip(*np.where(hit)):
        I, J = ia[i], ib[j]
        out.append((pP[I] + s_[i, j] * (pP[I + 1] - pP[I]), pQ[J] + t_[i, j] * (pQ[J + 1] - pQ[J])))
    return out


def refine(d, q, coeffs_fun, p0, b0):
    from scipy.optimize import fsolve

    def F(x):
        L = left_value(d, q, x[0])
        R = right_value(d, q, coeffs_fun, x[1])
        if L is None or R is None:
            return [1e3, 1e3]
        return L - R
    x, info, ier, msg = fsolve(F, [p0, b0], full_output=True, xtol=1e-13)
    return x, ier == 1, np.max(np.abs(F(x)))


def singular_amplitude(d, q):
    """psi_*(rho) = A rho^sigma (static singular solution w = A|x|^sigma), A^{q+1} = 1/(sigma(sigma+d-2))."""
    sig = 2.0 / (q + 1)
    if sig + d - 2 <= 0:
        return None
    return (1.0 / (sig * (sig + d - 2))) ** (1.0 / (q + 1))


def predicted_ratio(d, q):
    """Oscillation of the linearisation at psi_* near rho=0: rho^{-(d-2)/2 +- i nu},
    nu^2 = q sigma (sigma+d-2) - (d-2)^2/4.  Consecutive profiles: psi_n(0)/psi_{n+1}(0) -> exp(sigma pi/nu)."""
    sig = 2.0 / (q + 1)
    nu2 = q * sig * (sig + d - 2) - (d - 2) ** 2 / 4.0
    if nu2 <= 0:
        return None
    return np.exp(sig * np.pi / np.sqrt(nu2))


def profile_on_grid(d, q, p, b, rr):
    """Evaluate the matched smooth profile (left branch on [0,RHO_M], right branch on [RHO_M,1])."""
    sig = 2.0 / (q + 1)
    resonant, cf = right_series(d, q)
    rhs = rhs_factory(d, q)
    A = (p ** (-q) - sig * (1 - sig) * p) / (2 * d)
    r0 = 1e-6
    solL = solve_ivp(rhs, [r0, RHO_M], [p + A * r0 ** 2, 2 * A * r0], method="DOP853", rtol=1e-12,
                     atol=1e-14, dense_output=True)
    co = cf(b)
    eps = 2e-3
    Y0 = [sum(cj * eps ** j for j, cj in enumerate(co)), -sum(j * cj * eps ** (j - 1) for j, cj in enumerate(co) if j)]
    solR = solve_ivp(rhs, [1 - eps, RHO_M], Y0, method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
    out = np.empty_like(rr)
    for k, r in enumerate(rr):
        if r <= r0:
            out[k] = p + A * r * r
        elif r <= RHO_M:
            out[k] = solL.sol(r)[0]
        elif r < 1 - eps:
            out[k] = solR.sol(r)[0]
        else:
            y = 1 - r
            out[k] = sum(cj * y ** j for j, cj in enumerate(co))
    return out


def find_profiles(d, q, pmin):
    sig = 2.0 / (q + 1)
    c = (1 / (sig * (1 - sig))) ** (1.0 / (q + 1))
    resonant, cf = right_series(d, q)
    ps = c * np.exp(np.linspace(np.log(pmin), np.log(50.0), 2000))
    if resonant:
        bs = np.linspace(-60.0, 60.0, 3001) * c
    else:
        bs = c * np.exp(np.linspace(np.log(0.02), np.log(50.0), 1500))
    PL = curve(lambda p: left_value(d, q, p), ps)
    PR = curve(lambda b: right_value(d, q, cf, b), bs)
    X = seg_intersections(PL, PR, ps, bs)
    sols = []
    for (p0, b0) in X:
        x, conv, res = refine(d, q, cf, p0, b0)
        if conv and res < 1e-8 and not any(abs(x[0] - s_[0]) < 1e-6 * x[0] for s_ in sols):
            sols.append(x)
    sols.sort(key=lambda x: -x[0])
    return c, resonant, sols, PL, PR


def main():
    cases = [(d, q) for q in (2, 3, 4) for d in (1, 2, 3)]
    results = {}
    for (d, q) in cases:
        results[(d, q)] = find_profiles(d, q, 0.02 if d == 1 else 1e-3)
    # ---------------- report
    for (d, q) in cases:
        c, resonant, sols, _, _ = results[(d, q)]
        A = singular_amplitude(d, q)
        ratio = predicted_ratio(d, q)
        print(f"d={d} q={q} c={c:.10f} resonant={resonant}: {len(sols)} smooth positive profile(s) "
              f"with psi(0) >= {0.02 if d == 1 else 1e-3}*c")
        if A is not None:
            print(f"    singular profile A rho^sigma: A/c = {A / c:.6f};  predicted psi_n(0)/psi_(n+1)(0) -> {ratio:.4f}")
        prev = None
        for x in sols:
            extra = "" if prev is None else f"   ratio to previous: {prev / x[0]:.4f}"
            lab = "psi(1)/c" if not resonant else "psi'(1)"
            val = x[1] / c if not resonant else x[1]
            print(f"    psi(0)/c = {x[0] / c:.10f},  {lab} = {val:.10f}{extra}")
            prev = x[0]
    # ---------------- figure
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9))
    d, q = 3, 2
    c, resonant, sols, PL, PR = results[(d, q)]
    ax = axes[0]
    ax.plot(PL[:, 0] / c, PL[:, 1] / c, '-', color='#1f5fa8', lw=1.0, label=r'regular at $\rho=0$')
    ax.plot(PR[:, 0] / c, PR[:, 1] / c, '-', color='#c4531b', lw=1.0, label=r'regular at $\rho=1$')
    for x in sols:
        v = left_value(d, q, x[0])
        ax.plot(v[0] / c, v[1] / c, 'ko', ms=4)
    ax.set_xlim(0.2, 1.4)
    ax.set_ylim(-0.6, 1.2)
    ax.set_xlabel(r'$\psi(1/2)/c$')
    ax.set_ylabel(r"$\psi'(1/2)/c$")
    ax.set_title(r'(a) shooting, $d=3$, $q=2$', fontsize=10)
    ax.legend(fontsize=7, loc='upper right')
    rr = np.linspace(0, 1, 801)
    for ax, d in zip(axes[1:], (2, 3)):
        c, resonant, sols, _, _ = results[(d, q)]
        for k, x in enumerate(sols[:4]):
            prof = profile_on_grid(d, q, x[0], x[1], rr)
            ax.plot(rr, prof / c, lw=1.4, label=(r'$\psi_{%d}$' % k) + (' (ODE)' if k == 0 else ''))
        A = singular_amplitude(d, q)
        ax.plot(rr, A * rr ** (2.0 / (q + 1)) / c, 'k--', lw=1.0, label=r'$A\rho^{\sigma}$ (singular)')
        ax.set_xlabel(r'$\rho=|x|/(T-t)$')
        ax.set_ylabel(r'$\psi(\rho)/c$')
        ax.set_ylim(0, 1.1)
        ax.set_title(r'(%s) smooth profiles, $d=%d$, $q=2$' % ('b' if d == 2 else 'c', d), fontsize=10)
        ax.legend(fontsize=7, loc='lower right')
    fig.tight_layout()
    fig.savefig(FIG)
    print("figure written to", os.path.normpath(FIG))


if __name__ == "__main__":
    main()
