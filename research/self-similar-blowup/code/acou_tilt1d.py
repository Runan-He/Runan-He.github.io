#!/usr/bin/env python3
"""1D numerics for model (b):  u_tt - u_xx = (u_t)^2,  u(0) = 0,  u_t(0) = phi(x).

Purpose: test the 'tilted ODE blowup' picture of Section 7.
  (i)  the blowup curve t = T(x) is spacelike (|T'| < 1) near its minimum;
  (ii) pointwise rate:  u_t(t,x) ~ c(x)/(T(x)-t)  with  c(x) = 1 - T'(x)^2
       (NOT c = 1 unless T'(x) = 0);
  (iii) self-similar profile at a point x1 of the curve with a = T'(x1):
       rho u_t(T(x1)-rho, x1+rho xi) -> (1-a^2)/(1+a xi),   rho -> 0,
       with an O(rho) error (spectral gap 1).

Scheme: exact characteristic-diamond identity for u_tt - u_xx = f with dt = dx = h,
  U^{n+1}_j + U^{n-1}_j - U^n_{j+1} - U^n_{j-1} = h^2 f(t_n, x_j),
f = u_t^2 with u_t = (U^{n+1}_j - U^{n-1}_j)/(2h) (implicit, solved exactly:
D = U^{n+1}-U^{n-1} = 2(1 - sqrt(1-K)), K = U^n_{j+1}+U^n_{j-1}-2U^{n-1}_j).
Second order; respects the light cone exactly (no boundary conditions: the computational
interval shrinks by one node per step on each side).  If 1-K < 0 the node is declared
blown up (NaN); NaNs then fill exactly the complement of the maximal development.

Run: python3 acou_tilt1d.py [--fig]     (writes figures/acou_tilt.pdf with --fig)
"""
import sys
import os
import numpy as np


def phi(x):
    return 2.0/(1.0 + x**2)


def run(h, L, tmax, store_every_x=20, slab=None):
    """Return (t_grid, x_sub, ut_sub[n, :]) and optional full-resolution slabs.
    slab = (x_lo, x_hi): store u_t at full resolution for x in [x_lo, x_hi] at all steps."""
    x = np.arange(-L, L + h/2, h)
    J = x.size
    u0 = np.zeros(J)
    u1 = phi(x)
    # third-order start: u_tt(0) = u0'' + u1^2, u_ttt(0) = u1'' + 2 u1 u_tt(0)
    u1xx = np.zeros(J)
    u1xx[1:-1] = (u1[2:] - 2*u1[1:-1] + u1[:-2])/h**2
    utt0 = u1**2
    uttt0 = u1xx + 2*u1*utt0
    Um = u0.copy()
    U = u0 + h*u1 + h**2/2*utt0 + h**3/6*uttt0
    U[0] = U[-1] = np.nan
    nmax = int(tmax/h)
    sub = np.arange(0, J, store_every_x)
    ts, uts = [], []
    if slab is not None:
        js = np.where((x >= slab[0]) & (x <= slab[1]))[0]
        slabs = []
    lo, hi = 1, J - 1                      # valid interior range for U (current level)
    for n in range(1, nmax):
        K = np.full(J, np.nan)
        K[lo:hi] = U[lo+1:hi+1] + U[lo-1:hi-1] - 2*Um[lo:hi]
        disc = 1.0 - K
        with np.errstate(invalid="ignore"):
            D = np.where(disc >= 0, 2.0*(1.0 - np.sqrt(np.abs(disc))), np.nan)
        Up = Um + D
        Up[:lo] = np.nan
        Up[hi:] = np.nan
        ut = D/(2*h)                       # u_t at level n (time n h)
        ts.append(n*h)
        uts.append(ut[sub].copy())
        if slab is not None:
            slabs.append(ut[js].copy())
        Um, U = U, Up
        lo, hi = lo + 1, hi - 1
        if np.all(np.isnan(ut[sub][(x[sub] > -1.5) & (x[sub] < 1.5)])):
            break
    out = (np.array(ts), x[sub], np.array(uts))
    if slab is not None:
        out = out + (x[js], np.array(slabs))
    return out


def fit_blowup(ts, ut, rho_window=(0.01, 0.04)):
    """For each x: fit 1/u_t = (T - t)/c on the last part of the time series.
    Two-step: rough T from the last finite time, then a fit on rho in rho_window."""
    nx = ut.shape[1]
    T = np.full(nx, np.nan)
    c = np.full(nx, np.nan)
    for i in range(nx):
        y = ut[:, i]
        good = np.isfinite(y) & (y > 0)
        if good.sum() < 50:
            continue
        Tr = ts[good][-1]          # the scheme breaks down at T - O(h)
        for _ in range(4):
            rho = Tr - ts
            sel = good & (rho > rho_window[0]) & (rho < rho_window[1])
            if sel.sum() < 10:
                break
            # include an O(rho^2) correction: 1/u_t = (T-t)/c + e (T-t)^2
            A = np.vstack([ts[sel], np.ones(sel.sum()), ts[sel]**2]).T
            coef = np.linalg.lstsq(A, 1/y[sel], rcond=None)[0]
            roots = np.roots([coef[2], coef[0], coef[1]])
            roots = roots[np.isreal(roots)].real
            if roots.size == 0:
                break
            Tr = roots[np.argmin(np.abs(roots - Tr))]
            slope = coef[0] + 2*coef[2]*Tr        # d(1/u_t)/dt at t = T
            c[i] = -1/slope
        T[i] = Tr
    return T, c


def main(make_fig=False):
    h = 1.0e-4
    L, tmax = 2.6, 1.2
    ts, xs, ut = run(h, L, tmax, store_every_x=50)
    T, c = fit_blowup(ts, ut)
    m = np.isfinite(T) & np.isfinite(c) & (np.abs(xs) < 1.2)
    xs_, T_, c_ = xs[m], T[m], c[m]
    Tp = np.gradient(T_, xs_)
    print("x      T(x)      T'(x)     c(x)     1-T'^2    diff")
    for i in range(0, xs_.size, max(1, xs_.size//16)):
        print(f"{xs_[i]:+.3f}  {T_[i]:.6f}  {Tp[i]:+.5f}  {c_[i]:.5f}  {1-Tp[i]**2:.5f}  {c_[i]-(1-Tp[i]**2):+.1e}")
    err = np.max(np.abs(c_ - (1 - Tp**2))[2:-2])
    print(f"max |c - (1-T'^2)| over |x|<1.2: {err:.2e};  min T at x = {xs_[np.argmin(T_)]:+.3f}, "
          f"T_min = {T_.min():.6f};  max |T'| = {np.max(np.abs(Tp)):.3f}")
    # profile at x1
    x1 = 0.8
    i1 = np.argmin(np.abs(xs_ - x1))
    x1 = xs_[i1]
    T1, a = T_[i1], Tp[i1]
    # second run with a full-resolution slab around x1
    ts2, _, _, xslab, uslab = run(h, L, tmax, store_every_x=500, slab=(x1 - 0.25, x1 + 0.25))
    rhos = [0.2, 0.1, 0.05, 0.025]
    profs = []
    for rho in rhos:
        n = np.argmin(np.abs(ts2 - (T1 - rho)))
        rho_n = T1 - ts2[n]
        xi = (xslab - x1)/rho_n
        sel = np.abs(xi) <= 1
        prof = rho_n*uslab[n, sel]
        target = (1 - a**2)/(1 + a*xi[sel])
        e = np.nanmax(np.abs(prof - target))
        e0 = np.nanmax(np.abs(prof - 1.0))
        profs.append((rho_n, xi[sel], prof))
        print(f"x1={x1:+.3f}, a=T'(x1)={a:+.4f}, rho={rho_n:.4f}: "
              f"sup|rho u_t - (1-a^2)/(1+a xi)| = {e:.2e},  sup|rho u_t - 1| = {e0:.2e}")
    if make_fig:
        figure(xs_, T_, Tp, c_, x1, a, profs)


def figure(xs, T, Tp, c, x1, a, profs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    ax[0].plot(xs, c, color="#2b6cb0", lw=1.6, label=r"fitted $c(x)$ in $u_t\sim c/(T(x)-t)$")
    ax[0].plot(xs, 1 - Tp**2, color="k", lw=0.8, ls="--", label=r"$1-T'(x)^2$")
    ax[0].plot(xs, T - T.min() + 0.0, color="#c05621", lw=1.0, label=r"$T(x)-T_{\min}$")
    ax[0].set_xlabel("$x$")
    ax[0].legend(frameon=False, fontsize=7, loc="lower left")
    ax[0].set_title("blowup curve and pointwise rate", fontsize=8.5)
    shades = ["#bfd3ea", "#8fb2da", "#4f86c6", "#1f4e8c"]
    for (rho, xi, pr), col in zip(profs, shades):
        ax[1].plot(xi, pr, color=col, lw=1.4, label=fr"$\rho={rho:.3f}$")
    xi = np.linspace(-1, 1, 200)
    ax[1].plot(xi, (1 - a**2)/(1 + a*xi), "k--", lw=0.9, label=r"$(1-a^2)/(1+a\xi)$")
    ax[1].axhline(1.0, color="0.6", lw=0.6, ls=":")
    ax[1].set_xlabel(r"$\xi=(x-x_1)/\rho$")
    ax[1].set_ylabel(r"$\rho\,\partial_t u$")
    ax[1].set_title(fr"profile at $x_1={x1:.2f}$, $a=T'(x_1)={a:.3f}$", fontsize=8.5)
    ax[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "..", "figures", "acou_tilt.pdf")
    fig.savefig(path)
    print("wrote", path)


if __name__ == "__main__":
    main("--fig" in sys.argv)
