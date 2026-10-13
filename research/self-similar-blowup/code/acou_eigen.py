#!/usr/bin/env python3
"""Numerical spectra for Section 7 (Chebyshev collocation, no hypergeometric input).

(b)  u_tt - Lap u = (u_t)^2, linearisation at the ODE blowup in similarity variables:
        lam g1 = -r g1' + g2
        lam g2 = Lap_l g1 - r g2' + g2,      Lap_l g = g'' + (d-1)/r g' - l(l+d-2)/r^2 g
     predicted:  lam = 1 - j,  j = l, l+1, l+2, ...   (j = degree of a polynomial free wave)
(W)  inviscid Westervelt degenerate blowup, perturbation of Q^2, rescaled ball |zeta| <= 1:
        lam g1 = -(3/4) r g1' + g2
        lam g2 = (9/16) Lap_l g1 - (3/4) r g2' + g2
     predicted:  lam = 1 - (3/4)(l+2n)  and  lam = -(3/4)(l+2n).

Discretisation: write g_i = r^l h_i(r) Y_l (so Lap_l g = r^l (h'' + (2l+d-1)/r h'),
Lam g = r^l (l h + r h')), Chebyshev points on [-1,1] with N odd (r = 0 is not a node),
h even (Trefethen, Spectral Methods in MATLAB, Ch. 11), unknowns at the nodes r_j > 0
including r = 1.  No boundary condition at r = 1 (the light cone / sound horizon
is characteristic and the operator degenerates there).

Numerical caveat (observed): the collocation matrices are extremely non-normal.  In double
precision the eigenvalues are accurate (~1e-8) only for N <~ 31; for N >~ 61 rounding errors
produce spurious complex pairs creeping in from Re(lam) ~ -2.  In high precision (mpmath,
40-60 digits) the eigenvalues at N = 41, 61 agree with the exact values to > 20 digits,
which shows that these pairs are rounding artefacts (pseudospectra), not eigenvalues.

Run: python3 acou_eigen.py          (double precision table, N = 21, 25, 31)
     python3 acou_eigen.py --mp     (+ high-precision check, N = 41 and 61; a few minutes)
     python3 acou_eigen.py --fig    (+ writes figures/acou_spectrum.pdf from the mp data)
"""
import sys
import os
import numpy as np
import mpmath as mp


# ------------------------------------------------------------------ double precision
def cheb(N):
    x = np.cos(np.pi*np.arange(N + 1)/N)
    c = np.hstack([2, np.ones(N - 1), 2])*(-1)**np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1/c)/(dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


COEF = {"b": (1.0, 1.0), "W": (0.75, 9/16)}   # (transport speed, Laplacian coefficient)


def operator(N, d, l, model="b"):
    """Matrix of the linearised operator on (g1, g2) at the positive nodes."""
    assert N % 2 == 1
    D, x = cheb(N)
    D2 = D @ D
    M = (N + 1)//2
    r = x[:M]

    def fold(A):                 # h(-r) = h(r)
        return A[:M, :M] + A[:M, M:][:, ::-1]
    Dr, D2r = fold(D), fold(D2)
    R, Ri, I = np.diag(r), np.diag(1/r), np.eye(M)
    Lap = D2r + (2*l + d - 1)*Ri @ Dr
    a, b = COEF[model]
    Lam = l*I + R @ Dr
    return np.block([[-a*Lam, I], [b*Lap, -a*Lam + I]])


def spectrum(N, d, l, model="b"):
    ev = np.linalg.eigvals(operator(N, d, l, model))
    return ev[np.argsort(-ev.real)]


# ------------------------------------------------------------------ high precision
def operator_mp(N, d, l, model="b"):
    x = [mp.cos(mp.pi*j/N) for j in range(N + 1)]
    c = [(2 if j in (0, N) else 1)*(-1)**j for j in range(N + 1)]
    D = mp.zeros(N + 1, N + 1)
    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                D[i, j] = mp.mpf(c[i])/c[j]/(x[i] - x[j])
    for i in range(N + 1):
        D[i, i] = -mp.fsum(D[i, j] for j in range(N + 1) if j != i)
    D2 = D*D
    M = (N + 1)//2
    a, b = (mp.mpf(1), mp.mpf(1)) if model == "b" else (mp.mpf(3)/4, mp.mpf(9)/16)
    A = mp.zeros(2*M, 2*M)
    for i in range(M):
        r = x[i]
        for j in range(M):
            Dr = D[i, j] + D[i, N - j]
            D2r = D2[i, j] + D2[i, N - j]
            lap = D2r + (2*l + d - 1)/r*Dr
            lam_ = r*Dr + (l if i == j else 0)
            A[i, j] = -a*lam_
            A[M + i, j] = b*lap
            A[M + i, M + j] = -a*lam_ + (1 if i == j else 0)
            if i == j:
                A[i, M + j] = 1
    return A


def spectrum_mp(N, d, l, model="b", dps=50):
    with mp.workdps(dps):
        E = mp.eig(operator_mp(N, d, l, model), left=False, right=False)
        E = sorted(E, key=lambda z: -mp.re(z))
        return np.array([complex(z) for z in E]), E


def predicted(d, l, model="b", lmin=-3.0):
    out = []
    if model == "b":
        j = l
        while 1 - j > lmin:
            out.append(1.0 - j)
            j += 1
    else:
        n = 0
        while True:
            a = 1 - 0.75*(l + 2*n)
            b = -0.75*(l + 2*n)
            if a <= lmin and b <= lmin:
                break
            out += [v for v in (a, b) if v > lmin]
            n += 1
        out = sorted(out, reverse=True)
    return np.array(out)


CASES = [("b", 1, 0), ("b", 1, 1), ("b", 3, 0), ("b", 3, 1), ("b", 3, 2),
         ("W", 1, 0), ("W", 1, 1), ("W", 3, 0), ("W", 3, 1), ("W", 3, 2)]


def report_double(lmin=-3.2):
    print("=== double precision, N = 17, 21, 25: eigenvalues with Re > %.1f ===" % lmin)
    for model, d, l in CASES:
        specs = [spectrum(N, d, l, model) for N in (17, 21, 25)]
        ref = specs[-1][specs[-1].real > lmin]
        pred = predicted(d, l, model, lmin)
        err = max(np.min(np.abs(pred - z)) for z in ref) if len(ref) else np.nan
        spread = max(max(np.min(np.abs(s - z)) for s in specs[:-1]) for z in ref) if len(ref) else np.nan
        print(f"{model} d={d} l={l}: N=25 {np.round(ref.real, 6)} | exact {pred} | "
              f"max err {err:.1e}, N-spread {spread:.1e}, max|Im| {np.max(np.abs(ref.imag)):.1e}")


def report_mp(lmin=-3.2, Ns=((41, 40), (61, 60))):
    print("=== high precision (mpmath); errors computed in working precision ===")
    data = {}
    for model, d, l in CASES:
        for N, dps in Ns:
            ev, E = spectrum_mp(N, d, l, model, dps)
            pred = predicted(d, l, model, lmin)
            with mp.workdps(dps):
                sel = [z for z in E if mp.re(z) > lmin]
                ok = len(sel) == len(pred)
                if ok:
                    # exact values as mp numbers (pred are quarter-integers -> exact)
                    exact = [mp.mpf(str(v)) for v in pred]
                    err = max(abs(z - e) for z, e in zip(sel, exact))
                    errs = mp.nstr(err, 3)
                else:
                    errs = "n/a"
            print(f"{model} d={d} l={l} N={N} dps={dps}: {np.round(ev[ev.real > lmin].real, 12)}  "
                  f"exact {pred}  count {'=' if ok else '!='}  max err {errs}")
            data[(model, d, l, N)] = ev
    return data


def figure(path, data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
    col = {0: "#2b6cb0", 1: "#c05621", 2: "#2f855a"}
    lmin = -4.2
    for ax, model, title in zip(axes, ("b", "W"),
                                (r"model (b): $u_{tt}-\Delta u=(u_t)^2$",
                                 r"inviscid Westervelt (degenerate blowup)")):
        rows = [(1, 0), (1, 1), (3, 0), (3, 1), (3, 2)]
        for y, (d, l) in enumerate(rows):
            # double precision N = 81 (rounding artefacts) in light grey
            sd = spectrum(81, d, l, model)
            sd = sd[(sd.real > lmin) & (np.abs(sd.imag) < 8)]
            ax.scatter(sd.real, np.full(sd.shape, y + 0.22), s=5, color="0.7", lw=0, zorder=1)
            ev = data[(model, d, l, 61)]
            ev = ev[ev.real > lmin]
            ax.scatter(ev.real, np.full(ev.shape, y), s=14, color=col[l], lw=0, zorder=3)
            pr = predicted(d, l, model, lmin)
            ax.scatter(pr, np.full(pr.shape, y), s=46, facecolors="none", edgecolors="k",
                       lw=0.7, zorder=4)
        ax.axvline(0, color="0.5", lw=0.6, zorder=0)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([f"$d={d},\\ \\ell={l}$" for d, l in rows])
        ax.set_xlabel(r"$\lambda$")
        ax.set_xlim(lmin, 1.4)
        ax.set_ylim(-0.6, len(rows) - 0.3)
        ax.set_title(title, fontsize=8.5)
    fig.text(0.5, 0.005, "filled: collocation $N=61$, 60 digits;  circles: exact;  "
             "grey (offset): double precision $N=81$ (rounding artefacts of a non-normal matrix)",
             ha="center", fontsize=7, color="0.25")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path)
    print("wrote", path)


if __name__ == "__main__":
    report_double()
    if "--mp" in sys.argv or "--fig" in sys.argv:
        data = report_mp()
        if "--fig" in sys.argv:
            here = os.path.dirname(os.path.abspath(__file__))
            figure(os.path.join(here, "..", "figures", "acou_spectrum.pdf"), data)
