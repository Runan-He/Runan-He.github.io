#!/usr/bin/env python3
"""
num_spectra.py -- independent spectral computation for the linearisation of
u_tt - Delta u = F(u) (or (u_t)^2) around the ODE blowup profile, in the
similarity coordinates  tau = -log(T-t), xi = x/(T-t)  of NOTATION.md.

First-order system (eigenmode e^{lam tau}(phi1,phi2)), power-type scaling s:
    lam phi1 = -Lambda phi1 - s phi1 + phi2
    lam phi2 = Delta phi1 - Lambda phi2 - (s+1) phi2 + V0 phi1
derivative model u_tt - Delta u = u_t^2, psi1 = u + log(T-t), psi2 = (T-t)u_t:
    lam phi1 = -Lambda phi1 + phi2
    lam phi2 = Delta phi1 - Lambda phi2 + phi2
Spherical-harmonic sector l:  phi_i = rho^l Y_l h_i(z), z = rho^2 in [0,1]:
    Lambda -> l + 2 z d/dz,   Delta -> 4 z d^2/dz^2 + (4l + 2d) d/dz.
Chebyshev collocation in z on N+1 Gauss-Lobatto points, NO boundary condition
at z = 1 (light cone, characteristic) and none at z = 0 (regularity is built
into the ansatz rho^l h(rho^2) with h polynomial).

Outputs (run from project root):  python3 code/num_spectra.py
  * stdout: tables of converged eigenvalues (Re lam >= -3) for all cases,
    comparison with lam = a_i - s - l - 2n, exact (rational) multiplicities,
    eigenfunction checks against symmetry generators and 2F1 polynomials.
  * figures/num_spectra.pdf: double vs high precision spectra (spurious branch).
Runtime: about 3-5 minutes on 4 cores (mpmath part is the bottleneck).
"""
import os
import sys
import time
from fractions import Fraction as Fr
from multiprocessing import Pool

import numpy as np
import mpmath as mp
import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
DPS = 30           # working precision of the high-precision runs
NMP = (12, 20)     # high-precision resolutions (exact for degree <= N)
NDP = (20, 40, 80)  # double-precision resolutions
REMIN = -3.0       # report eigenvalues with Re >= REMIN


# --------------------------------------------------------------------------
# model data:  (name, params, s, V0, deriv, dims)
# --------------------------------------------------------------------------
def cases():
    out = []
    for q in (2, 3, 4):
        sig = Fr(2, q + 1)
        out.append(("MEMS", "q=%d" % q, -sig, q * sig * (1 - sig), False, (1, 2, 3)))
    out.append(("exp", "", Fr(0), Fr(2), False, (1, 2, 3)))
    for p in (3, 5, 7):
        s = Fr(2, p - 1)
        out.append(("power", "p=%d" % p, s, p * s * (s + 1), False, (3,)))
    out.append(("deriv", "", Fr(0), Fr(0), True, (1, 3)))
    return out


def sectors(d):
    return (0, 1) if d == 1 else (0, 1, 2)


def formula(s, V0, deriv, l, remin=REMIN):
    """closed-form guess of NOTATION.md (for deriv: the corrected version)."""
    if deriv:
        roots = [Fr(1), Fr(0)]           # lam = 1-l-2n, -l-2n
        shift = Fr(0)
    else:
        disc = 1 + 4 * V0
        r = sp.sqrt(sp.Rational(disc.numerator, disc.denominator))
        roots = [(-1 + r) / 2, (-1 - r) / 2]
        roots = [Fr(str(sp.nsimplify(x))) for x in roots]
        shift = s
    lams = []
    for a in roots:
        n = 0
        while True:
            lam = a - shift - l - 2 * n
            if lam < remin - Fr(1, 10 ** 9):
                break
            lams.append(lam)
            n += 1
    return sorted(lams, reverse=True)


# --------------------------------------------------------------------------
# Chebyshev matrices (double and mpmath)
# --------------------------------------------------------------------------
def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def opmat(N, d, l, s, V0, deriv):
    D, x = cheb(N)
    z = (1 + x) / 2
    Dz = 2 * D
    Z = np.diag(z)
    I = np.eye(N + 1)
    Lam = l * I + 2 * Z @ Dz
    Lap = 4 * Z @ Dz @ Dz + (4 * l + 2 * d) * Dz
    if deriv:
        A = np.block([[-Lam, I], [Lap, -Lam + I]])
    else:
        A = np.block([[-Lam - s * I, I], [Lap + V0 * I, -Lam - (s + 1) * I]])
    return A, z


def opmat_mp(N, d, l, s, V0, deriv):
    x = [mp.cos(mp.pi * j / N) for j in range(N + 1)]
    c = [(2 if j in (0, N) else 1) * (-1) ** j for j in range(N + 1)]
    n = N + 1
    D = mp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            if i != j:
                D[i, j] = mp.mpf(c[i]) / c[j] / (x[i] - x[j])
    for i in range(n):
        D[i, i] = -mp.fsum(D[i, j] for j in range(n) if j != i)
    z = [(1 + xx) / 2 for xx in x]
    Dz = 2 * D
    D2 = Dz * Dz
    s = mp.mpf(s.numerator) / s.denominator
    V0 = mp.mpf(V0.numerator) / V0.denominator
    A = mp.zeros(2 * n, 2 * n)
    for i in range(n):
        for j in range(n):
            lam = 2 * z[i] * Dz[i, j]
            A[i, j] = -lam
            A[n + i, n + j] = -lam
            A[n + i, j] = 4 * z[i] * D2[i, j] + (4 * l + 2 * d) * Dz[i, j]
        if deriv:
            A[i, i] += -l
            A[n + i, n + i] += -l + 1
        else:
            A[i, i] += -l - s
            A[n + i, i] += V0
            A[n + i, n + i] += -l - s - 1
        A[i, n + i] += 1
    return A, z


def mp_eigs(args):
    N, d, l, s, V0, deriv = args
    mp.mp.dps = DPS
    A, _ = opmat_mp(N, d, l, s, V0, deriv)
    ev = mp.eig(A, left=False, right=False)
    return [complex(e) for e in ev], [mp.nstr(mp.re(e), 28) for e in ev]


# --------------------------------------------------------------------------
# exact rational computation on the invariant polynomial subspace P_n x P_n
# (monomial basis); gives exact eigenvalues and Jordan structure.
# --------------------------------------------------------------------------
def exact_jordan(d, l, s, V0, deriv, nmax=5):
    n = nmax + 1
    S = sp.Rational(s.numerator, s.denominator)
    V = sp.Rational(V0.numerator, V0.denominator)
    A = sp.zeros(2 * n, 2 * n)
    for k in range(n):
        lamk = l + 2 * k
        if deriv:
            A[k, k] = -lamk
            A[n + k, n + k] = -lamk + 1
        else:
            A[k, k] = -lamk - S
            A[n + k, n + k] = -lamk - S - 1
            A[n + k, k] = V
        A[k, n + k] = 1
        if k >= 1:  # Laplacian: z^k -> k(4k-4+4l+2d) z^{k-1}
            A[n + k - 1, k] = k * (4 * k - 4 + 4 * l + 2 * d)
    res = {}
    for lam, alg in A.eigenvals().items():
        if lam < REMIN - sp.Rational(1, 10 ** 9):
            continue
        M = A - lam * sp.eye(2 * n)
        geo = 2 * n - M.rank()
        res[lam] = (alg, geo)
    return res


# --------------------------------------------------------------------------
def hyp_poly(a, b, c, z):
    """2F1(a,b;c;z) for a in -N0 (terminating)."""
    tot, term, k = 0.0, 1.0, 0
    while True:
        tot += term
        if a + k == 0:
            break
        term *= (a + k) * (b + k) / ((c + k) * (k + 1)) * z
        k += 1
        if k > 50:
            break
    return tot


def eigvec_checks():
    """eigenvectors (double precision, N = 20) vs symmetry generators / 2F1."""
    lines = []
    N = 20
    for (name, par, s, V0, deriv, dims) in cases():
        for d in dims:
            for l in sectors(d):
                A, z = opmat(N, d, l, float(s), float(V0), deriv)
                w, V = np.linalg.eig(A)
                fl_all = formula(s, V0, deriv, l, remin=-1.5)
                for lamf in sorted(set(fl_all), reverse=True):
                    i = np.argmin(abs(w - float(lamf)))
                    h1 = V[: N + 1, i]
                    h2 = V[N + 1:, i]
                    k0 = np.argmax(abs(h1))
                    h1, h2 = h1 / h1[k0], h2 / h1[k0]
                    # expected h1: 2F1 polynomial in z (alpha or beta in -N0)
                    if deriv:
                        a = float(lamf) - 1.0
                        a1, a2 = 0.0, -1.0
                    else:
                        a = float(lamf + s)
                        r = np.sqrt(1 + 4 * float(V0))
                        a1, a2 = (-1 + r) / 2, (-1 - r) / 2
                    al, be = (a + l - a1) / 2, (a + l - a2) / 2
                    if abs(al - round(al)) < 1e-9 and round(al) <= 0:
                        hex_ = np.array([hyp_poly(round(al), be, d / 2 + l, zz) for zz in z])
                    else:
                        hex_ = np.array([hyp_poly(round(be), al, d / 2 + l, zz) for zz in z])
                    hex_ = hex_ / hex_[k0]
                    err1 = np.max(abs(h1 - hex_))
                    # phi2 = (lam + s + Lambda) phi1 (resp. lam + Lambda for deriv)
                    D, _ = cheb(N)
                    Lh = l * h1 + 2 * z * (2 * D @ h1)
                    h2ex = (float(lamf) + (0 if deriv else float(s))) * h1 + Lh
                    err2 = np.max(abs(h2 - h2ex))
                    lines.append((name, par, d, l, lamf, float(w[i].real), err1, err2,
                                  h2[0].real / h1[0].real, fl_all.count(lamf) > 1))
    return lines


# --------------------------------------------------------------------------
def fmt_fr(x):
    x = Fr(x)
    return str(x.numerator) if x.denominator == 1 else "%d/%d" % (x.numerator, x.denominator)


def main():
    t0 = time.time()
    allcases = cases()
    jobs = []
    for (name, par, s, V0, deriv, dims) in allcases:
        for d in dims:
            for l in sectors(d):
                for N in NMP:
                    jobs.append((N, d, l, s, V0, deriv))
    print("high-precision runs: %d matrices, dps=%d" % (len(jobs), DPS), flush=True)
    with Pool(min(4, os.cpu_count() or 1)) as pool:
        mpres = pool.map(mp_eigs, jobs, chunksize=1)
    mpdict = {j[:3] + (j[3], j[5]): r for j, r in zip(jobs, mpres)}
    print("  done in %.0f s" % (time.time() - t0), flush=True)

    print("\n=== converged eigenvalues (Re >= %g) ===" % REMIN)
    print("high precision (dps=%d): Chebyshev N=%d and N=%d agree to < 1e-8" % ((DPS,) + NMP))
    print("dbl80: number of these reproduced to 1e-6 by double precision at N=20/40/80")
    summary = []
    for (name, par, s, V0, deriv, dims) in allcases:
        for d in dims:
            for l in sectors(d):
                evA = np.array(mpdict[(NMP[0], d, l, s, deriv)][0])
                evB = np.array(mpdict[(NMP[1], d, l, s, deriv)][0])
                conv = []
                for e in evA:
                    if e.real < REMIN - 1e-6:
                        continue
                    # converged in N: N=12 and N=20 agree (1e-8 tolerates the
                    # O(eps^(1/2)) splitting of a 2x2 Jordan block at dps=30)
                    if np.min(abs(evB - e)) < 1e-8:
                        conv.append(e)
                conv = sorted(conv, key=lambda e: -e.real)
                fl = formula(s, V0, deriv, l)
                # match with multiplicity
                fl_f = [float(x) for x in fl]
                err = 0.0
                unmatched = list(conv)
                for lamf in fl_f:
                    if not unmatched:
                        err = np.inf
                        break
                    k = int(np.argmin([abs(u - lamf) for u in unmatched]))
                    err = max(err, abs(unmatched[k] - lamf))
                    unmatched.pop(k)
                ex = exact_jordan(d, l, s, V0, deriv)
                jord = ["%s(alg %d,geo %d)" % (fmt_fr(Fr(str(k))), a, g)
                        for k, (a, g) in sorted(ex.items(), key=lambda t: -t[0]) if a != g]
                # double precision
                dbl = []
                for N in NDP:
                    A, _ = opmat(N, d, l, float(s), float(V0), deriv)
                    ev = np.linalg.eigvals(A)
                    good = sum(np.min(abs(ev - f)) < 1e-6 for f in set(fl_f))
                    dbl.append(good)
                summary.append((name, par, d, l, conv, fl, err, unmatched, jord, dbl,
                                len(set(fl_f))))
                print("%-6s %-4s d=%d l=%d : %s" % (name, par, d, l,
                      ", ".join("%.12g" % e.real if abs(e.imag) < 1e-12 else
                                "%.6g%+.2gi" % (e.real, e.imag) for e in conv)))
                print("        formula: %s | max|num-formula| = %.1e | extra: %s"
                      % (", ".join(fmt_fr(x) for x in fl), err,
                         [complex(np.round(u, 8)) for u in unmatched]))
                print("        exact (rational) Jordan blocks: %s | double-prec. resolved"
                      " (of %d distinct) N=20/40/80: %s" % (jord or "none (all semisimple)",
                                                            len(set(fl_f)), dbl))

    # Jordan splitting seen in high precision
    print("\n=== splitting of double eigenvalues in high precision (N=%d) ===" % NMP[1])
    for (name, par, d, l, conv, fl, err, unm, jord, dbl, nd) in summary:
        if jord:
            for lamf in sorted(set(fl), reverse=True):
                if fl.count(lamf) == 2:
                    cl = sorted([e for e in conv if abs(e - float(lamf)) < 1e-6],
                                key=lambda e: e.real)
                    print("  %s %s d=%d l=%d lam=%s : |split| = %.2e, |mean-lam| = %.1e"
                          " (dps=%d: split ~ 10^-15 = sqrt(eps) signals a Jordan block)"
                          % (name, par, d, l, fmt_fr(lamf),
                             abs(cl[-1] - cl[0]) if len(cl) == 2 else np.nan,
                             abs(np.mean(cl) - float(lamf)) if cl else np.nan, DPS))

    print("\n=== eigenfunction checks (double precision, N=20) ===")
    print("err1 = max|h1 - 2F1 polynomial|, err2 = max|h2 - (lam+s+Lambda)h1|,"
          " ratio = h2/h1 at z=1")
    checks = eigvec_checks()
    for (name, par, d, l, lamf, lamn, e1, e2, ratio, dbl) in checks:
        if lamf in (1, 0):
            tag = ""
            if lamf == 1 and l == 0:
                tag = "time translation"
            if lamf == 0 and l == 1:
                tag = "Lorentz boost (phi1 ~ xi)"
            if lamf == 0 and l == 0:
                tag = "shift u -> u + const"
            print("  %-6s %-4s d=%d l=%d lam=%-3s num=%+.3e err1=%.1e err2=%.1e h2/h1=%.6f  %s"
                  % (name, par, d, l, fmt_fr(lamf), lamn - float(lamf), e1, e2, ratio, tag))
    ok = [max(c[6], c[7]) for c in checks if not c[9]]
    print("  all simple eigenvalues with lam >= -1.5: max eigenvector deviation from the"
          " 2F1 polynomial = %.1e (%d checks)" % (max(ok), len(ok)))
    for c in checks:
        if c[9]:
            print("  (double eigenvalue %s %s d=%d l=%d lam=%s: double-precision eigenvector"
                  " deviation %.1e, as expected for a Jordan block)" % (c[0], c[1], c[2], c[3],
                                                                      fmt_fr(c[4]), max(c[6], c[7])))

    make_figure()
    print("\ntotal time %.0f s" % (time.time() - t0))


def make_figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10})
    panels = [("MEMS q=2, d=1, $\\ell=0$", 1, 0, Fr(-2, 3), Fr(4, 9), False),
              ("MEMS q=2, d=3, $\\ell=0$", 3, 0, Fr(-2, 3), Fr(4, 9), False),
              ("power p=3, d=3, $\\ell=0$", 3, 0, Fr(1), Fr(6), False)]
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True)
    cols = {20: "#1f77b4", 40: "#ff7f0e", 80: "#2ca02c"}
    mks = {20: "o", 40: "s", 80: "^"}
    with Pool(3) as pool:
        mpres = pool.map(mp_eigs, [(32, d, l, s, V0, der) for (_, d, l, s, V0, der) in panels])
    for ax, (title, d, l, s, V0, der), mpr in zip(axs, panels, mpres):
        for N in NDP:
            A, _ = opmat(N, d, l, float(s), float(V0), der)
            ev = np.linalg.eigvals(A)
            ax.plot(ev.real, ev.imag, mks[N], ms=4, mfc="none", color=cols[N],
                    label="double, N=%d" % N)
        evm = np.array(mpr[0])
        ax.plot(evm.real, evm.imag, "k.", ms=4, label="%d digits, N=32" % DPS)
        for lamf in formula(s, V0, der, l, remin=-8):
            ax.axvline(float(lamf), color="0.85", lw=0.6, zorder=0)
        lc = (d - 1) / 2 - float(s)
        ax.axvline(lc, color="r", ls="--", lw=0.8)
        ax.text(lc + 0.1, 5.2, r"$\frac{d-1}{2}-s$", color="r", fontsize=9)
        ax.set_xlim(-8, 2.2)
        ax.set_ylim(-6, 6)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(r"Re $\lambda$")
    axs[0].set_ylabel(r"Im $\lambda$")
    axs[0].legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "num_spectra.pdf"))
    print("wrote figures/num_spectra.pdf")


if __name__ == "__main__":
    main()
