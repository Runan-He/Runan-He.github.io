"""geom_wavemaps_spectrum.py -- eigenvalues of the linearisation around the
self-similar co-rotational wave map f0(rho) = 2 arctan(rho/sqrt(d-2)).

Run:  python3 geom_wavemaps_spectrum.py   (30-60 minutes, mostly the
argument-principle step; set QUICK=1 in the environment to skip it).
Writes figures/geom_spectrum.pdf (relative to the project root).

Three independent methods:
 (1) Chebyshev collocation in y = rho^2 (companion linearisation of the
     quadratic pencil), N = 40..120, to locate candidates and test convergence;
 (2) Frobenius series at y=0 and y=1 + Wronskian at y=1/2, normalised by
     1/Gamma(1-sigma2) so that it is entire in lambda; Newton refinement to
     25+ digits;
 (3) argument principle for the Wronskian on rectangles, to show that no
     eigenvalue in the box was missed.
Also: the SUSY partner problem (gauge mode removed) and the Sturm count of
bound states of the hyperbolic Schroedinger operator (no eigenvalues with
Re lambda > (d-1)/2).
"""
import os
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
import geom_modes as G

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures", "geom_spectrum.pdf")
QUICK = os.environ.get("QUICK", "0") == "1"


def collocation_table(d, Ns=(40, 60, 80, 100, 120), cut=-2.6):
    print(f"\n--- collocation, d = {d} (eigenvalues with Re > {cut}) ---")
    out = {}
    for N in Ns:
        ev = G.collocation_spectrum(d, N)
        ev = ev[ev.real > cut]
        out[N] = ev
        print(f"N={N:4d}: " + ", ".join(f"{z.real:+.10f}{z.imag:+.2e}j" for z in ev[:8]))
    return out


def refined(d, which="orig", guesses=(1.001, -0.55, -2.05)):
    vals = []
    for g in guesses:
        r = G.refine(d, g, which, dps=35)
        if all(abs(r - v) > 1e-10 for v in vals):
            vals.append(r)
    return vals


def sturm_nodes(d, eta_max=40.0):
    """Zero-energy (kappa = 0, i.e. lambda = (d-1)/2) solution of the 1D form
    -psi'' + Q psi = 0,  Q = (n-1)(n-3)/(4 sinh^2 eta) + V(tanh eta)/cosh^2 eta,
    regular at eta = 0 (psi ~ eta^{(n-1)/2}).  By Sturm oscillation theory its
    number of zeros on (0, inf) equals the number of bound states of
    H = -d^2/deta^2 + Q, i.e. of eigenvalues lambda > (d-1)/2.  Returns that number
    and psi(eta_max) (sign/size of the tail)."""
    n = d + 2
    m = (n - 1) / 2

    def Q(e):
        r = np.tanh(e)
        return (m * (m - 1)) / np.sinh(e) ** 2 - 8 * (d - 1) * (d - 2) / (d - 2 + r * r) ** 2 / np.cosh(e) ** 2

    def rhs(e, Y):
        return [Y[1], Q(e) * Y[0]]

    e0 = 1e-3
    sol = solve_ivp(rhs, [e0, eta_max], [e0 ** m, m * e0 ** (m - 1)], rtol=1e-11,
                    atol=1e-30, dense_output=True)
    es = np.linspace(e0, eta_max, 400001)
    psi = sol.sol(es)[0]
    return int(np.sum(np.sign(psi[1:]) != np.sign(psi[:-1]))), psi[-1]


def main():
    mp.mp.dps = 30
    results = {}
    for d in (3, 4, 5):
        tab = collocation_table(d)
        o = refined(d, "orig")
        s = refined(d, "susy", guesses=(1.001, -0.55, -2.05))
        results[d] = dict(tab=tab, orig=o, susy=s)
        print(f"refined eigenvalues (orig), d={d}: " + ", ".join(mp.nstr(z.real, 20) for z in o))
        print(f"refined eigenvalues (SUSY partner), d={d}: " + ", ".join(mp.nstr(z.real, 20) for z in s))
        nz, tail = sturm_nodes(d)
        print(f"Sturm: zero-energy solution has {nz} zeros on (0,inf)  => "
              f"{nz} eigenvalues with Re lambda > (d-1)/2 = {(d-1)/2}")

    # first two stable eigenvalues as a function of d
    print("\n--- lambda_1(d), lambda_2(d) (refined) ---")
    dd = list(range(3, 13))
    lam1, lam2 = [], []
    g1, g2 = -0.55, -2.05
    for d in dd:
        a = G.refine(d, g1, "orig", dps=30).real
        b = G.refine(d, g2, "orig", dps=30).real
        lam1.append(float(a))
        lam2.append(float(b))
        g1, g2 = a, b
        print(f"d={d:2d}: lambda_1 = {mp.nstr(a, 15)}, lambda_2 = {mp.nstr(b, 15)}")

    if not QUICK:
        print("\n--- argument principle (number of zeros of the entire Wronskian) ---")
        for d in (3, 4, 5):
            for which in ("orig", "susy"):
                box = (-1.9, 4.0, -15.0, 15.0)
                cnt, mx = G.count_zeros(d, box, which, npts=1200)
                print(f"d={d} {which:4s} box Re in [-1.9,4], Im in [-15,15]: "
                      f"winding = {cnt:.4f} (max phase step {mx:.2f} rad)")
            # strip just left of the big box (contains lambda_2 only)
            cnt, mx = G.count_zeros(d, (-2.25, -1.95, -6.0, 6.0), "orig", npts=500)
            print(f"d={d} orig box Re in [-2.25,-1.95], Im in [-6,6]: "
                  f"winding = {cnt:.4f} (max phase step {mx:.2f} rad)")

    # ------------------------------------------------------------ figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "font.family": "serif"})
    cols = {3: "#2a78d6", 4: "#eb6834", 5: "#1baf7a"}
    ink, muted = "#0b0b0b", "#52514e"
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw={"width_ratios": [1.25, 1]})
    a = ax[0]
    a.axvspan(-2.0, 1.6, color="#e9e8e4", lw=0, zorder=0)
    a.axvline(0, color=muted, lw=0.6, ls=":")
    ev60 = G.collocation_spectrum(3, 60)
    ev120 = G.collocation_spectrum(3, 120)
    a.plot(ev60.real, ev60.imag, "o", ms=4.5, mfc="none", mec=cols[3], mew=0.9,
           label="collocation $N=60$", zorder=3)
    a.plot(ev120.real, ev120.imag, "x", ms=4, color=muted, mew=0.9,
           label="collocation $N=120$", zorder=3)
    tru = [float(z.real) for z in results[3]["orig"]]
    a.plot(tru, [0] * len(tru), "o", ms=7, color=cols[3], mec="white", mew=1.2,
           label="Wronskian (exact to 25 digits)", zorder=4)
    for x, lab in zip(tru, ["$1$ (gauge)", "$-0.54247$", "$-2$"]):
        a.annotate(lab, (x, 0), xytext=(0, 9), textcoords="offset points", ha="center",
                   fontsize=8, color=ink)
    a.set_xlim(-4.2, 1.6)
    a.set_ylim(-6, 6)
    a.set_xlabel(r"Re $\lambda$")
    a.set_ylabel(r"Im $\lambda$")
    a.set_title(r"(a) $d=3$: spectrum in the $\lambda$-plane", fontsize=9, loc="left")
    a.text(-1.95, -5.4, "Re$\\,\\lambda>-2$", fontsize=8, color=muted)
    a.legend(loc="upper right", fontsize=7, frameon=False)
    b = ax[1]
    b.axhline(1, color=muted, lw=0.8, ls="--")
    b.text(12, 1.08, r"gauge mode $\lambda=1$", ha="right", fontsize=8, color=muted)
    b.plot(dd, lam1, "-o", color=cols[3], ms=4.5, lw=2, label=r"$\lambda_1(d)$")
    b.plot(dd, lam2, "-s", color=cols[4], ms=4.5, lw=2, label=r"$\lambda_2(d)$")
    b.set_xlabel(r"space dimension $d$")
    b.set_ylabel(r"$\lambda$")
    b.set_ylim(-2.9, 1.4)
    b.set_xticks(dd)
    b.legend(loc="center left", frameon=False, fontsize=8, bbox_to_anchor=(0.02, 0.42))
    b.set_title(r"(b) first stable eigenvalues vs. $d$", fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig(FIG)
    print("wrote", os.path.normpath(FIG))


if __name__ == "__main__":
    main()
