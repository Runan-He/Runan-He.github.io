"""geom_membrane.py -- relativistic membranes (timelike extremal hypersurfaces,
Nambu-Goto / scalar Born-Infeld): symbolic and numerical checks for Section 6.

 1. Graph equation in R^{1+(d+1)}: z = u(t,x),  L = sqrt(1 - u_t^2 + |grad u|^2);
    non-divergence form; scaling u -> mu^{-1} u(mu t, mu x)  (s = -1).
 2. Radial self-similar ansatz u = (T-t) U(rho), rho = r/(T-t): every solution
    analytic at rho = 0 is constant (moving flat plane) -- Taylor recursion.
 3. The "lightlike self-similar solutions" u = +-(T-t) sqrt(1-rho^2) are null
    cones: they solve the non-divergence form but L = 0 (not timelike).
 4. Spherical k-membrane R(t,omega) in R^{k+1}: L = -R^k sqrt(1 - R_t^2 + |grad_S R|^2/R^2)
    (induced-metric determinant checked for k = 2); energy E = R^k/sqrt(1-R_t^2);
    k = 1 exact solution R = E cos(t/E); collapse asymptotics
    R = s - s^{2k+1}/(2(2k+1)E^2) + ..., s = T - t.
 5. Linearisation about the collapsing sphere, delta = sum delta_l(t) Y_l:
      delta'' - 2k (R0'/R0) delta' + (mu_l - k) R0^{2k-2} E^{-2} delta = 0,
      mu_l = l(l+k-1);  exponents at s = 0 are 0 and 2k+1, so generically
      delta_l(T) != 0 and delta_l / R0 ~ (T-t)^{-1}: relative deformations grow for
      every l (lambda = 1 for every l in similarity variables).
    Numerical integration for k = 2, l = 0..4.
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

t, T, r, rho, s = sp.symbols("t T r rho s", positive=True)


def graph_equation(d):
    u = sp.Function("u")(t, r)
    L = sp.sqrt(1 - sp.diff(u, t) ** 2 + sp.diff(u, r) ** 2)
    EL = sp.diff(r ** (d - 1) * sp.diff(u, t) / L, t) - sp.diff(r ** (d - 1) * sp.diff(u, r) / L, r)
    ut, ur = sp.diff(u, t), sp.diff(u, r)
    utt, utr, urr = sp.diff(u, t, 2), sp.diff(u, t, r), sp.diff(u, r, 2)
    # claimed non-divergence form (radial):
    ND = ((1 + ur**2) * utt - 2 * ut * ur * utr - (1 - ut**2) * urr
          - (1 - ut**2 + ur**2) * (d - 1) / r * ur)
    chk = sp.simplify(EL * L**3 / r ** (d - 1) - ND)
    return u, ND, chk


def radial_selfsimilar(d, order=8):
    u, ND, chk = graph_equation(d)
    print(f"1. d={d}: divergence form * L^3 == non-divergence form:", chk)
    U = sp.Function("U")
    ND_ss = ND.subs(u, (T - t) * U(r / (T - t))).doit()
    ND_ss = sp.simplify((ND_ss * (T - t)).subs(r, rho * (T - t)))
    # Taylor recursion around rho = 0 with U(0) = c, U even
    c = sp.symbols("c")
    a = sp.symbols("a1:%d" % (order + 1))
    Useries = c + sum(a[j] * rho ** (2 * (j + 1)) for j in range(order // 2))
    ex = ND_ss.subs(U(rho), Useries).doit()
    ex = sp.series(sp.simplify(ex), rho, 0, order).removeO()
    sol = sp.solve([sp.expand(ex).coeff(rho, k) for k in range(0, order - 2, 2)], a[: order // 2 - 1], dict=True)
    print(f"2. d={d}: Taylor coefficients of regular radial self-similar solutions:", sol)
    # 3. null cone
    for sgn in (1, -1):
        un = sgn * sp.sqrt((T - t) ** 2 - r**2)
        val = sp.simplify(ND.subs(u, un).doit())
        Lsq = sp.simplify(1 - sp.diff(un, t) ** 2 + sp.diff(un, r) ** 2)
        print(f"3. d={d}: null cone sign {sgn}: residual {val}, 1-u_t^2+u_r^2 = {Lsq}")


def sphere_checks():
    th, ph = sp.symbols("theta phi")
    R = sp.Function("R")(t, th, ph)
    om = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
    X = R * om
    Xt, Xa, Xb = X.diff(t), X.diff(th), X.diff(ph)
    vecs = [Xt, Xa, Xb]
    eta = sp.diag(-1, 1, 1, 1)
    g = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            g[i, j] = -1 * 0 + (vecs[i].T * vecs[j])[0] - (1 if (i == 0 and j == 0) else 0)
    detg = sp.simplify(g.det())
    grad2 = sp.diff(R, th) ** 2 + sp.diff(R, ph) ** 2 / sp.sin(th) ** 2
    claim = -(R**4) * sp.sin(th) ** 2 * (1 - sp.diff(R, t) ** 2 + grad2 / R**2)
    print("4. k=2 induced metric: det g - claim =", sp.simplify(detg - claim))
    E = sp.symbols("E", positive=True)
    Rk1 = E * sp.cos(t / E)
    print("4. k=1 exact solution: R/sqrt(1-R_t^2) - E =",
          sp.simplify(Rk1 / sp.sqrt(1 - sp.diff(Rk1, t) ** 2) - E).subs(t, E / 3))
    for k in (1, 2, 3):
        c1 = sp.Rational(1, 2 * (2 * k + 1))
        Rs = s - c1 * s ** (2 * k + 1) / E**2
        Rt = -sp.diff(Rs, s)
        en = sp.series(Rs**k / sp.sqrt(1 - Rt**2), s, 0, 2 * k + 1).removeO()
        print(f"4. k={k}: energy of asymptotic expansion = {sp.simplify(en)} (should be E + O(s^(2k)))")


def linearised_sphere_symbolic():
    k, E, mu = sp.symbols("k E mu", positive=True)
    eps, Y, gY2 = sp.symbols("epsilon Y gY2")
    R0 = sp.Function("R0")(t)
    dl = sp.Function("delta")(t)
    R = R0 + eps * dl * Y
    A = 1 - sp.diff(R, t) ** 2 + eps**2 * dl**2 * gY2 / R**2
    L2 = sp.series(-(R**k) * sp.sqrt(A), eps, 0, 3).removeO().coeff(eps, 2)
    L2 = sp.expand(L2).subs(gY2, mu).subs(Y**2, 1)
    EL = sp.diff(L2, dl) - sp.diff(sp.diff(L2, sp.diff(dl, t)), t)
    R0d = -sp.sqrt(1 - R0 ** (2 * k) / E**2)
    R0dd = -k * R0 ** (2 * k - 1) / E**2
    EL = EL.subs(sp.Derivative(R0, (t, 2)), R0dd).subs(sp.Derivative(R0, t), R0d)
    claim = sp.diff(dl, t, 2) - 2 * k * (R0d / R0) * sp.diff(dl, t) + (mu - k) * R0 ** (2 * k - 2) / E**2 * dl
    # replace functions by symbols and show EL = (nonvanishing factor) * claim
    a0, a1, a2, x = sp.symbols("a0 a1 a2 x", positive=True)
    rep = lambda e: e.subs(sp.Derivative(dl, (t, 2)), a2).subs(sp.Derivative(dl, t), a1).subs(dl, a0).subs(R0, x)
    factor = sp.simplify(rep(EL).coeff(a2) / rep(claim).coeff(a2))
    print("5. linearised equation: EL = F * claim with F =", factor,
          "; remainder:", sp.simplify(sp.expand(rep(EL) - factor * rep(claim))))


def linearised_sphere_numeric(k=2, E=1.0, ells=(0, 1, 2, 3, 4), data=(1.0, 0.0)):
    """Background released from rest at R0(0) = E^{1/k}; integrate background and
    delta_l with delta(0)=1, delta'(0)=0 (and for l=0 the time-translation mode)."""
    R_init = E ** (1.0 / k)

    def rhs(tt, Y, mu):
        R, Rd, d0, d1 = Y
        Rdd = -k * R ** (2 * k - 1) / E**2
        return [Rd, Rdd, d1, 2 * k * Rd / R * d1 - (mu - k) * R ** (2 * k - 2) / E**2 * d0]

    stop = lambda tt, Y, mu: Y[0] - 1e-5 * R_init
    stop.terminal = True
    out = {}
    for l in ells:
        mu = l * (l + k - 1)
        sol = solve_ivp(rhs, [0, 10], [R_init, 0, data[0], data[1]], args=(mu,), rtol=1e-12, atol=1e-15,
                        events=stop, dense_output=True)
        Tc = sol.t_events[0][0]
        ts = np.linspace(0, Tc, 4000)
        Y = sol.sol(ts)
        out[l] = (Tc, ts, Y)
        print(f"5. data={data}, k={k}, l={l}: collapse at t~{Tc:.6f}; delta_l(T) ~ {Y[2, -1]:+.6f}; "
              f"delta_l/R0 at R0=1e-5: {Y[2, -1] / Y[0, -1]:+.3e}")
    return out


if __name__ == "__main__":
    for d in (1, 2, 3):
        radial_selfsimilar(d)
    sphere_checks()
    linearised_sphere_symbolic()
    linearised_sphere_numeric(data=(1.0, 0.0))
    linearised_sphere_numeric(data=(0.0, 1.0))
    linearised_sphere_numeric(k=1, data=(0.0, 1.0))
