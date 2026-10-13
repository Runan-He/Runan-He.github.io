#!/usr/bin/env python3
"""Light numerical check of the spectrum of the linearisation around self-similar MEMS touchdown
profiles  w = (T-t)^sigma psi(x/(T-t)),  w_tt - Delta w = -lam w^{-q}  (lam = 1 by scaling).

(A) ODE profile psi = c (constant potential V0 = 2q(q-1)/(q+1)^2):
    eigenvalue ODE in the sector of spherical harmonics of degree l, g = r^l h(z), z = r^2:
       z(1-z) h'' + (l + d/2 - (b + 3/2) z) h' - (b^2 + b - V)/4 h = 0,   b = mu - sigma + l,
    a quadratic eigenvalue problem in mu, solved by Chebyshev collocation on z in [0,1]
    (no boundary conditions: the collocation space consists of functions smooth at z = 0, 1).
    Prediction (Section 3):  mu = 1 - l - 2n  and  mu = -(2 - 2 sigma) - l - 2n.
(B) Independent check for d = 1: Chebyshev collocation of the full first-order operator
       L(f1,f2) = (-xi f1' + sigma f1 + f2,  f1'' - xi f2' - (1-sigma) f2 + V0 f1)  on [-1,1].
(C) Excited radial profiles psi_n (d = 3, q = 2; found by code/mems_profiles_shooting.py):
    same collocation with V(r) = q psi_n(r)^{-q-1}; count eigenvalues with Re mu > 0 per sector.
Writes figures/mems_spectrum.pdf.
"""
import os
import sys
import numpy as np
import scipy.linalg as sla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
FIG = os.path.join(HERE, "..", "figures", "mems_spectrum.pdf")


def cheb(N):
    """Chebyshev differentiation matrix on x_j = cos(pi j/N) (Trefethen)."""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    cc = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(cc, 1 / cc) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(1))
    return D, x


def sector_eigs(d, q, l, N, Vfun=None):
    """Eigenvalues mu of the sector-l problem; Vfun(z) gives the potential (default: constant V0)."""
    sig = 2.0 / (q + 1)
    V0 = 2.0 * q * (q - 1) / (q + 1) ** 2
    D, x = cheb(N)
    z = (1 - x) / 2          # z in [0,1]
    Dz = -2 * D
    Dz2 = Dz @ Dz
    V = V0 * np.ones_like(z) if Vfun is None else Vfun(z)
    I = np.eye(N + 1)
    Z = np.diag(z)
    gam = l + d / 2.0
    A0 = np.diag(z * (1 - z)) @ Dz2 + np.diag(gam - 1.5 * z) @ Dz + np.diag(V / 4.0)
    A1 = -Z @ Dz - 0.25 * I
    A2 = -0.25 * I
    # A0 h + b A1 h + b^2 A2 h = 0 ; companion linearisation in b
    Zr = np.zeros_like(I)
    Abig = np.block([[Zr, I], [-A0, -A1]])
    Bbig = np.block([[I, Zr], [Zr, A2]])
    b = sla.eig(Abig, Bbig, right=False)
    b = b[np.isfinite(b)]
    return b - l + sig


def converged(ev1, ev2, tol=2e-3, window=-1.9):
    """Eigenvalues of ev1 in Re mu > window that are reproduced by ev2 (finer grid) up to tol;
    clusters closer than 1e-3 (split multiple eigenvalues) are merged.  Moderate N is used on
    purpose: for large N the collocation also resolves the non-smooth (H^k, small k) generalised
    eigenfunctions at z = 1 and spurious eigenvalues from Re mu <= -omega_k drift into the window."""
    out = []
    for e in ev1:
        if e.real < window:
            continue
        if np.min(np.abs(ev2 - e)) < tol * max(1, abs(e)):
            if not any(abs(e - o) < 1e-3 for o in out):
                out.append(e)
    return np.array(sorted(out, key=lambda v: -v.real))


def predicted(q, l, window=-1.9):
    sig = 2.0 / (q + 1)
    out = []
    for n in range(10):
        for m in (1 - l - 2 * n, -(2 - 2 * sig) - l - 2 * n):
            if m >= window:
                out.append(m)
    return np.array(sorted(set(np.round(out, 12)), reverse=True))


def full_system_d1(q, N):
    sig = 2.0 / (q + 1)
    V0 = 2.0 * q * (q - 1) / (q + 1) ** 2
    D, x = cheb(N)
    I = np.eye(N + 1)
    X = np.diag(x)
    L = np.block([[-X @ D + sig * I, I], [D @ D + V0 * I, -X @ D - (1 - sig) * I]])
    return np.linalg.eigvals(L)


def main():
    report = []
    # (A) ODE profile, sectors
    allok = True
    for q in (2, 3, 4):
        for d in (1, 2, 3):
            for l in ((0, 1) if d == 1 else (0, 1, 2)):
                e1 = sector_eigs(d, q, l, 20)
                e2 = sector_eigs(d, q, l, 28)
                ev = converged(e1, e2)
                pr = predicted(q, l)
                ok = (len(ev) == len(pr)) and np.allclose(np.sort(ev.real), np.sort(pr), atol=1e-4) \
                    and np.all(np.abs(ev.imag) < 1e-3)
                allok &= ok
                report.append(f"q={q} d={d} l={l}: computed {np.round(ev.real, 8)}  predicted {np.round(pr, 8)}  "
                              f"{'OK' if ok else 'MISMATCH'}")
    print("\n".join(report))
    print("(A) all sectors agree with the hypergeometric prediction:", allok)

    # (B) full first-order operator in d = 1
    for q in (2, 3, 4):
        sig = 2.0 / (q + 1)
        e1 = full_system_d1(q, 40)
        e2 = full_system_d1(q, 56)
        ev = converged(e1, e2)
        pr = sorted(set(np.round(list(predicted(q, 0)) + list(predicted(q, 1)), 10)), reverse=True)
        print(f"(B) d=1 full system q={q}: computed {np.round(ev.real, 7)}  (imag max {np.max(np.abs(ev.imag)) if len(ev) else 0:.1e})")
        print(f"    predicted (both parities)   {np.round(pr, 7)}")

    # (C) excited profiles, d = 3, q = 2.  psi_n varies on the scale psi_n(0)^{3/2} near rho = 0,
    #     so psi_2 needs N >= 120 (spurious complex pairs at N <= 80 are discretisation artefacts).
    import mems_profiles_shooting as ms
    d, q = 3, 2
    sig = 2.0 / (q + 1)
    c = (1 / (sig * (1 - sig))) ** (1.0 / (q + 1))
    resonant, cf = ms.right_series(d, q)
    guesses = [(0.2453363549, 0.5695749308, (40, 80), (0, 1, 2)), (0.0565466455, 0.5860270951, (120, 160), (0,))]
    excited = []
    rr = np.concatenate([np.linspace(0, 0.05, 20001)[:-1], np.linspace(0.05, 1, 8001)])
    for p0, b0, Ns, sectors in guesses:
        x, conv, res = ms.refine(d, q, cf, p0 * c, b0 * c)
        assert conv and res < 1e-9
        prof = ms.profile_on_grid(d, q, x[0], x[1], rr)
        Vfun = (lambda z, prof=prof: q * np.interp(np.sqrt(z), rr, prof) ** (-q - 1))
        spec = {}
        for l in sectors:
            e1 = sector_eigs(d, q, l, Ns[0], Vfun)
            e2 = sector_eigs(d, q, l, Ns[1], Vfun)
            ev = converged(e1, e2, tol=2e-3, window=-0.3)
            spec[l] = ev
            print(f"(C) excited profile psi(0)/c={x[0] / c:.4f}, l={l}, N={Ns}: eigenvalues with Re mu > -0.3: "
                  f"{np.round(ev, 4)}")
        excited.append((x, spec))

    # ---------------- figure: ODE profile spectrum (d=3, q=2) + excited profile spectra (l=0)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.9))
    ax = axes[0]
    d, q = 3, 2
    cols = {0: '#1f5fa8', 1: '#c4531b', 2: '#3b8c3b'}
    for l in (0, 1, 2):
        ev = converged(sector_eigs(d, q, l, 20), sector_eigs(d, q, l, 28))
        pr = predicted(q, l)
        ax.plot(pr, np.full_like(pr, l), 'o', mfc='none', mec='k', ms=10)
        ax.plot(ev.real, np.full_like(ev.real, l), 'x', color=cols[l], ms=8, mew=1.8,
                label=r'$\ell=%d$' % l)
    ax.axvline(0, color='0.6', lw=0.8)
    ax.axvline(-2.0 / 3, color='0.6', lw=0.8, ls=':')
    ax.set_yticks([0, 1, 2])
    ax.set_ylim(-0.6, 2.6)
    ax.set_xlim(-2.3, 1.3)
    ax.set_xlabel(r'eigenvalue $\mu$ (all real)')
    ax.set_ylabel(r'sector $\ell$')
    ax.set_title(r'(a) ODE profile, $d=3$, $q=2$: $\times$ collocation, $\circ$ formula', fontsize=10)
    ax.annotate(r'$-\omega_0=-2/3$', xy=(-2.0 / 3, 2.35), xytext=(-1.55, 2.35), fontsize=8,
                arrowprops=dict(arrowstyle='->', lw=0.6))
    ax.legend(fontsize=8, loc='lower left')
    ax = axes[1]
    ev0 = converged(sector_eigs(d, q, 0, 20), sector_eigs(d, q, 0, 28), window=-1.2)
    rows = [(r'$\psi_0=c$', ev0)] + [(r'$\psi_{%d}$' % (k + 1), sp[0]) for k, (x, sp) in enumerate(excited)]
    for k, (lab, ev) in enumerate(rows):
        ev = ev[ev.real > -1.2]
        ax.plot(ev.real, np.full(len(ev), k), 'o', color=['k', '#1f5fa8', '#c4531b'][k], ms=7)
        for e in ev:
            ax.annotate('%.3g' % e.real, (e.real, k), textcoords='offset points', xytext=(0, 7),
                        ha='center', fontsize=7)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.set_xscale('symlog', linthresh=1.0)
    ax.set_xlim(-1.3, 30)
    ax.set_xticks([-1, 0, 1, 3, 10, 20])
    ax.set_xticklabels(['-1', '0', '1', '3', '10', '20'])
    ax.axvline(0, color='0.6', lw=0.8)
    ax.set_xlabel(r'eigenvalue $\mu$ (symlog scale; all real)')
    ax.set_title(r'(b) radial sector $\ell=0$, $d=3$, $q=2$: ODE vs excited profiles', fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG)
    print("figure written to", os.path.normpath(FIG))


if __name__ == "__main__":
    main()
