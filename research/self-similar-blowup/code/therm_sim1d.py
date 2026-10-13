#!/usr/bin/env python3
"""therm_sim1d.py -- numerical experiments for Section 4 (hyperbolic thermal explosion).

(1) Spatially homogeneous ODEs in physical units (tau_r, delta):
      reduced   tau th'' + th' = delta e^th
      full      tau th'' + th' = delta e^th + tau delta e^th th'
    -> th + 2 log(T-t) -> log(2 tau/delta)  vs  th + log(T-t) -> -log(delta).

(2) 1D nonlinear evolution in similarity variables (Chebyshev collocation on
    xi in [-1,1], no boundary conditions -- the cone boundary is outflow), RK4 in
    tau, blowup time T found by bisection (time-translation instability):
      reduced  u_tt + g u_t - u_xx = e^u       psi1 = u + 2 log(T-t), psi2 = (T-t)u_t
      full     (d_t + g)(u_t - e^u) = u_xx     psi1 = u +   log(T-t), psi2 = (T-t)u_t
    (a) reduced, g = 0, data from an explicit Liouville solution: compare T, the
        boost parameter and the e^{-tau} rate with the exact solution;
    (b) full, g = 0 and g = 1, generic data: convergence to a tilted profile
        Psi_v = (log(1-v^2) - log(1+v xi), 1/(1+v xi)) at rate ~ tau e^{-tau}
        (Jordan block at lambda = -1), with v = slope of the blowup curve.

Outputs figures/therm_blowup.pdf and prints a summary.
Run time: a few minutes.
"""
import os
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', 'figures', 'therm_blowup.pdf')


# ---------------------------------------------------------------------------
# (1) homogeneous ODEs
# ---------------------------------------------------------------------------
def ode_blowup(model, tau_r=0.5, delta=1.0, th0=0.0, dth0=0.0, thmax=40.0):
    def f(t, Y):
        th, dth = Y
        src = delta * np.exp(th)
        if model == 'reduced':
            return [dth, (src - dth) / tau_r]
        return [dth, (src + tau_r * src * dth - dth) / tau_r]
    ev = lambda t, Y: Y[0] - thmax
    ev.terminal = True
    sol = solve_ivp(f, [0, 100], [th0, dth0], method='DOP853', rtol=1e-13, atol=1e-13,
                    events=ev, dense_output=True)
    tend = sol.t[-1]
    th_end, dth_end = sol.y[:, -1]
    # blowup time from the leading asymptotics th' ~ k/(T-t), k = 2 (reduced) or 1 (full)
    k = 2.0 if model == 'reduced' else 1.0
    T = tend + k / dth_end
    tt = np.linspace(0, tend, 4000)
    th = sol.sol(tt)[0]
    return T, tt, th


# ---------------------------------------------------------------------------
# (2) similarity-variable solver
# ---------------------------------------------------------------------------
def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def cc_weights(N):
    """Clenshaw-Curtis quadrature weights on the Chebyshev points (Trefethen, clencurt)."""
    theta = np.pi * np.arange(N + 1) / N
    w = np.zeros(N + 1)
    ii = np.arange(1, N)
    v = np.ones(N - 1)
    if N % 2 == 0:
        w[0] = w[N] = 1 / (N**2 - 1)
        for k in range(1, N // 2):
            v -= 2 * np.cos(2 * k * theta[ii]) / (4 * k**2 - 1)
        v -= np.cos(N * theta[ii]) / (N**2 - 1)
    else:
        w[0] = w[N] = 1 / N**2
        for k in range(1, (N - 1) // 2 + 1):
            v -= 2 * np.cos(2 * k * theta[ii]) / (4 * k**2 - 1)
    w[ii] = 2 * v / N
    return w


class Sim:
    def __init__(self, model, gamma=0.0, N=32, dt=None):
        self.model, self.gamma, self.N = model, gamma, N
        self.D, self.xi = cheb(N)
        self.D2 = self.D @ self.D
        self.w = cc_weights(N)
        self.dt = dt if dt is not None else 0.6 / N**2
        self.k = 2.0 if model == 'reduced' else 1.0   # psi1 = u + k log(T-t)

    def rhs(self, tau, p1, p2, T):
        D, D2, xi, g = self.D, self.D2, self.xi, self.gamma
        dp1 = -xi * (D @ p1) + p2 - self.k
        if self.model == 'reduced':
            dp2 = D2 @ p1 - xi * (D @ p2) - p2 + np.exp(p1) - g * T * np.exp(-tau) * p2
        else:
            e = np.exp(p1)
            dp2 = D2 @ p1 - xi * (D @ p2) - p2 + e * p2 - g * T * np.exp(-tau) * (p2 - e)
        return dp1, dp2

    def initial(self, u0, u1, T, x0=0.0):
        X = x0 + T * self.xi
        return u0(X) + self.k * np.log(T), T * u1(X)

    def run(self, u0, u1, T, x0=0.0, taumax=12.0, ref=None, store_every=None, blow=2.0):
        """Evolve; return ('up'|'down'|'ok', taus, snapshots). 'up' = psi1 grows (T too late)."""
        p1, p2 = self.initial(u0, u1, T, x0)
        if ref is None:
            ref = p1[self.N // 2] if self.N % 2 == 0 else p1.mean()
        ref0 = p1[self.N // 2]
        nsteps = int(round(taumax / self.dt))
        dt = self.dt
        taus, snaps = [], []
        tau = 0.0
        if store_every is None:
            store_every = max(1, int(round(0.05 / dt)))
        for n in range(nsteps):
            k1 = self.rhs(tau, p1, p2, T)
            k2 = self.rhs(tau + dt / 2, p1 + dt / 2 * k1[0], p2 + dt / 2 * k1[1], T)
            k3 = self.rhs(tau + dt / 2, p1 + dt / 2 * k2[0], p2 + dt / 2 * k2[1], T)
            k4 = self.rhs(tau + dt, p1 + dt * k3[0], p2 + dt * k3[1], T)
            p1 = p1 + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            p2 = p2 + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            tau += dt
            if n % store_every == 0:
                taus.append(tau)
                snaps.append((p1.copy(), p2.copy()))
            # psi2(0) -> 2 (reduced) or 1 (full) on every boosted/tilted profile; if T is too
            # early, (T-t)u_t -> 0; if T is too late, psi2 blows up at finite tau.
            dev = p2[self.N // 2] / (2.0 if self.model == 'reduced' else 1.0) - 1.0
            if not np.isfinite(dev) or dev > blow:
                return 'up', np.array(taus), snaps
            if dev < -0.7:
                return 'down', np.array(taus), snaps
        self.final_dev = p2[self.N // 2] / (2.0 if self.model == 'reduced' else 1.0) - 1.0
        return 'ok', np.array(taus), snaps

    def find_T(self, u0, u1, Tlo, Thi, x0=0.0, taumax=12.0, iters=48):
        """Bisection on the blowup time of the point x0 using the unstable mode."""
        for _ in range(iters):
            Tm = 0.5 * (Tlo + Thi)
            status, _, _ = self.run(u0, u1, Tm, x0, taumax=taumax)
            if status == 'ok':      # threshold not reached: decide by the sign of the final deviation
                status = 'up' if self.final_dev > 0 else 'down'
            if status == 'up':
                Thi = Tm
            else:
                Tlo = Tm
            if Thi - Tlo < 1e-15:
                break
        return 0.5 * (Tlo + Thi)

    def h1norm(self, f):
        return np.sqrt(self.w @ (f**2 + (self.D @ f)**2))


def profile(model, v, xi, delta=1.0):
    if model == 'reduced':
        A0 = 1 / np.sqrt(1 - v**2)
        return np.log(2 / delta) - 2 * np.log(A0 * (1 + v * xi)), 2 / (1 + v * xi)
    return np.log(1 - v**2) - np.log(1 + v * xi), 1 / (1 + v * xi)


# ---------------------------------------------------------------------------
# Liouville test data
# ---------------------------------------------------------------------------
def liouville_data(eps=0.08, T0=1.0):
    sig, eta, t, x = sp.symbols('sigma eta t x', real=True)
    F = (T0 - sig) / 2 + eps * (sp.sin(2 * sig) / 2 + sig**2 / 4)
    G = -eta / 2 + eps * (sp.cos(3 * eta) / 3 - eta**2 / 8)
    Fp, Gp = sp.diff(F, sig), sp.diff(G, eta)
    u = sp.log(8 * Fp.subs(sig, t + x) * Gp.subs(eta, t - x) / (F.subs(sig, t + x) + G.subs(eta, t - x))**2)
    ut = sp.diff(u, t)
    u0 = sp.lambdify(x, u.subs(t, 0), 'numpy')
    u1 = sp.lambdify(x, ut.subs(t, 0), 'numpy')
    uf = sp.lambdify((t, x), u, 'numpy')
    utf = sp.lambdify((t, x), ut, 'numpy')
    # blowup point above x = 0 : F(T) + G(T) = 0
    Ts = sp.nsolve((F + G.subs(eta, sig)).subs(sig, t), t, T0)
    p = float(Fp.subs(sig, Ts)); q = float(Gp.subs(eta, Ts))
    v = (q - p) / (p + q)
    # curvature of the blowup curve (implicit differentiation)
    F2 = float(sp.diff(F, sig, 2).subs(sig, Ts)); G2 = float(sp.diff(G, eta, 2).subs(eta, Ts))
    Tpp = -(F2 * (1 + v)**2 + G2 * (1 - v)**2) / (p + q)
    return dict(u0=u0, u1=u1, u=uf, ut=utf, T=float(Ts), v=v, Tpp=Tpp)


def main():
    out = {}
    # ---------------- (1) ODEs ----------------
    tau_r, delta = 0.5, 1.0
    Tr, ttr, thr = ode_blowup('reduced', tau_r, delta)
    Tf, ttf, thf = ode_blowup('full', tau_r, delta)
    print(f"(1) homogeneous ODEs (tau_r={tau_r}, delta={delta}, th(0)=th'(0)=0):")
    print(f"    reduced: T = {Tr:.10f}; th+2log(T-t) at end = {thr[-1] + 2*np.log(Tr-ttr[-1]):.6f}"
          f"  (predicted log(2 tau/delta) = {np.log(2*tau_r/delta):.6f})")
    print(f"    full   : T = {Tf:.10f}; th+ log(T-t) at end = {thf[-1] + np.log(Tf-ttf[-1]):.6f}"
          f"  (predicted -log(delta) = {-np.log(delta):.6f})")
    # full model: exact reduction th' - delta e^th = C e^{-t/tau}, C = th'(0) - delta e^{th(0)} = -1
    # => y = e^{-th}: y' = -delta - C e^{-t/tau} y, y(0)=1 ; T is the zero of y
    from scipy.integrate import quad
    from scipy.optimize import brentq
    Cc = -1.0
    def y_of_t(t):
        sol = solve_ivp(lambda s, y: [-delta - Cc * np.exp(-s / tau_r) * y[0]], [0, t], [1.0],
                        rtol=1e-13, atol=1e-15)
        return sol.y[0, -1]
    Tf_exact = brentq(y_of_t, 0.5, 5.0, xtol=1e-14)
    print(f"    full   : T from the first-order reduction = {Tf_exact:.10f}")
    out['ode'] = (Tr, ttr, thr, Tf, ttf, thf, tau_r, delta)

    # ---------------- (2a) reduced, Liouville ----------------
    L = liouville_data()
    sim = Sim('reduced', gamma=0.0, N=32)
    Tnum = sim.find_T(L['u0'], L['u1'], L['T'] - 0.05, L['T'] + 0.05, taumax=14.0)
    print(f"(2a) reduced/Liouville: exact T = {L['T']:.12f}, bisection T = {Tnum:.12f}, "
          f"|diff| = {abs(Tnum-L['T']):.1e}; exact slope v = T'(0) = {L['v']:.6f}, T''(0) = {L['Tpp']:.4f}")
    status, taus, snaps = sim.run(L['u0'], L['u1'], L['T'], taumax=14.0)
    P1, P2 = profile('reduced', L['v'], sim.xi)
    err_num = np.array([sim.h1norm(a - P1) + sim.h1norm(b - P2) for a, b in snaps])
    # exact error along the same tau values
    err_ex = []
    for tau in taus:
        s = L['T'] * np.exp(-tau)
        xx = s * sim.xi
        a = L['u'](L['T'] - s, xx) + 2 * np.log(s)
        b = s * L['ut'](L['T'] - s, xx)
        err_ex.append(sim.h1norm(a - P1) + sim.h1norm(b - P2))
    err_ex = np.array(err_ex)
    # predicted leading term: s * g1 with g1 = T'' A0^2 (1-xi^2)/(1+v xi) for psi1
    sel = (taus > 3) & (taus < 9)
    slope_ex = np.polyfit(taus[sel], np.log(err_ex[sel]), 1)[0]
    slope_num = np.polyfit(taus[sel], np.log(err_num[sel]), 1)[0]
    print(f"     fitted decay rates on tau in [3,9]: exact {slope_ex:.4f}, numerical {slope_num:.4f} (predicted -1)")
    out['liou'] = (taus, err_ex, err_num)

    # boost parameter from the numerical profile vs the slope of the blowup curve
    p1e, p2e = snaps[np.argmin(abs(taus - 9.0))]
    v_fit = -(sim.D @ p2e)[sim.N // 2] / 2.0
    print(f"     boost from numerical profile at tau=9: v = {v_fit:.6f} (exact {L['v']:.6f})")

    # ---------------- (2b) full model, generic data ----------------
    T0 = 1.0
    u0 = lambda x: -np.log(T0) + 0.3 * np.exp(-2 * (x - 0.3)**2) + 0.1 * x
    u1 = lambda x: 1 / T0 + 0.2 * np.cos(2 * x) + 0.15 * np.sin(x)
    res_full = {}
    for g in [0.0, 1.0]:
        simF = Sim('full', gamma=g, N=32)
        Tst = simF.find_T(u0, u1, 0.3, 1.5, taumax=13.0)
        # slope of blowup curve by bisection at x0 = +-h
        h = 0.02
        Tp = simF.find_T(u0, u1, Tst - 0.1, Tst + 0.1, x0=h, taumax=11.0)
        Tm = simF.find_T(u0, u1, Tst - 0.1, Tst + 0.1, x0=-h, taumax=11.0)
        v_curve = (Tp - Tm) / (2 * h)
        Tpp_curve = (Tp - 2 * Tst + Tm) / h**2
        status, tausF, snF = simF.run(u0, u1, Tst, taumax=13.0)
        p1e, p2e = snF[np.argmin(abs(tausF - 11.0))]
        v_fit = -(simF.D @ p2e)[simF.N // 2]
        P1, P2 = profile('full', v_fit, simF.xi)
        errF = np.array([simF.h1norm(a - P1) + simF.h1norm(b - P2) for a, b in snF])
        print(f"(2b) full, gamma={g}: T = {Tst:.10f}; tilt from profile v = {v_fit:.6f}, "
              f"slope of blowup curve T'(0) = {v_curve:.6f}, T''(0) ~ {Tpp_curve:.3f}")
        sel = (tausF > 3) & (tausF < 9)
        sl = np.polyfit(tausF[sel], np.log(errF[sel]), 1)[0]
        sl2 = np.polyfit(tausF[sel], np.log(errF[sel] / tausF[sel]), 1)[0]
        print(f"     fit log(err) ~ {sl:.4f} tau ;  fit log(err/tau) ~ {sl2:.4f} tau  (Jordan block: err ~ tau e^-tau)")
        res_full[g] = (tausF, errF, v_fit, v_curve)
    out['full'] = res_full

    # reduced with damping (gamma = 1), same data as Liouville
    simD = Sim('reduced', gamma=1.0, N=32)
    TD = simD.find_T(L['u0'], L['u1'], 0.5, 2.0, taumax=13.0)
    status, tausD, snD = simD.run(L['u0'], L['u1'], TD, taumax=13.0)
    p1e, p2e = snD[np.argmin(abs(tausD - 11.0))]
    vD = -(simD.D @ p2e)[simD.N // 2] / 2.0
    P1, P2 = profile('reduced', vD, simD.xi)
    errD = np.array([simD.h1norm(a - P1) + simD.h1norm(b - P2) for a, b in snD])
    sel = (tausD > 3) & (tausD < 9)
    slD = np.polyfit(tausD[sel], np.log(errD[sel]), 1)[0]
    print(f"(2c) reduced, gamma=1 (Liouville data): T = {TD:.10f}, boost v = {vD:.5f}, decay rate {slD:.4f} (predicted -1)")
    out['damped'] = (tausD, errD)

    return out


def make_figure(out):
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9, 'legend.fontsize': 7.5})
    fig, ax = plt.subplots(1, 3, figsize=(10.0, 3.1))
    # panel 1: ODE asymptotics -- (T-t) theta_t -> 2 (reduced) and 1 (full)
    Tr, ttr, thr, Tf, ttf, thf, tau_r, delta = out['ode']
    for (T_, tt_, th_, col, lab) in [(Tr, ttr, thr, 'C0', 'reduced'), (Tf, ttf, thf, 'C3', 'full')]:
        dth = np.gradient(th_, tt_)
        m = (T_ - tt_) > 2e-5
        ax[0].semilogx((T_ - tt_)[m], ((T_ - tt_) * dth)[m], color=col, label=lab)
    ax[0].axhline(2, color='C0', ls=':', lw=0.8)
    ax[0].axhline(1, color='C3', ls=':', lw=0.8)
    ax[0].invert_xaxis()
    ax[0].set_ylim(0, 2.6)
    ax[0].set_xlabel(r'$T-t$')
    ax[0].set_ylabel(r'$(T-t)\,\theta_t$')
    ax[0].set_title(r'(a) homogeneous blowup, $\tau_r=%g,\ \delta=%g$' % (tau_r, delta))
    ax[0].legend(loc='center left')
    # panel 2: reduced (Liouville) rates
    taus, err_ex, err_num = out['liou']
    tausD, errD = out['damped']
    ax[1].semilogy(taus, err_ex, color='k', lw=1.6, label='exact (Liouville)')
    ax[1].semilogy(taus, err_num, color='C1', ls='--', label='numerical, $\\gamma=0$')
    ax[1].semilogy(tausD, errD, color='C2', ls='-.', label='numerical, $\\gamma=1$')
    ax[1].semilogy(taus, 0.5 * err_ex[np.argmin(abs(taus - 4))] * np.exp(-(taus - 4)), color='gray',
                   ls=':', label=r'$\propto e^{-\tau}$')
    ax[1].set_xlabel(r'$\tau=-\log(T-t)+\log T$')
    ax[1].set_title(r'(b) reduced: $\|\Psi(\tau)-\Psi_{a}\|$')
    ax[1].set_ylim(1e-7, 10)
    ax[1].legend(loc='upper right')
    # panel 3: e^{tau} * error: linear growth (Jordan block, full) vs. constant (reduced)
    for g, col in [(0.0, 'C3'), (1.0, 'C4')]:
        tausF, errF, v_fit, v_curve = out['full'][g]
        m = tausF < 12.5
        ax[2].plot(tausF[m], (np.exp(tausF) * errF)[m], color=col,
                   label=f'full, $\\gamma={g:g}$ ($v={v_fit:.3f}$)')
    m = taus < 12.5
    ax[2].plot(taus[m], (np.exp(taus) * err_ex)[m], color='k', label='reduced (Liouville)')
    m = tausD < 12.5
    ax[2].plot(tausD[m], (np.exp(tausD) * errD)[m], color='C2', ls='-.', label=r'reduced, $\gamma=1$')
    ax[2].set_xlabel(r'$\tau$')
    ax[2].set_title(r'(c) $e^{\tau}\,\|\Psi(\tau)-\Psi_{\rm lim}\|$')
    ax[2].legend(loc='upper left')
    for a in ax:
        a.grid(alpha=0.3, lw=0.5)
    fig.tight_layout()
    fig.savefig(FIG)
    print("wrote", os.path.normpath(FIG))


if __name__ == "__main__":
    import pickle
    cache = os.environ.get('THERM_CACHE')          # optional: reuse results for figure tweaks
    if cache and os.path.exists(cache):
        with open(cache, 'rb') as fh:
            out = pickle.load(fh)
    else:
        out = main()
        if cache:
            with open(cache, 'wb') as fh:
                pickle.dump(out, fh)
    make_figure(out)
