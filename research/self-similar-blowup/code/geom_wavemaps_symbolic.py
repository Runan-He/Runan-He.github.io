"""geom_wavemaps_symbolic.py -- symbolic checks for Section 6 (wave maps).

Checks (all must print 0 / True):
 1. f0 = 2 arctan(rho/sqrt(d-2)) solves the self-similar ODE, d = 3..10.
 2. u = r psi turns the co-rotational equation in R^{1+d} into a radial wave
    equation in R^{1+(d+2)} with a smooth nonlinearity.
 3. Mode equation from the PDE (w = (T-t)^{-1-lam} g(r/(T-t)), s = 1).
 4. Potential V(rho) = (d-1)(cos 2f0 - 1)/rho^2 = -8(d-1)(d-2)/(d-2+rho^2)^2.
 5. Gauge mode: lam = 1, g = f0'(rho) ~ 1/(d-2+rho^2).
 6. Hyperbolic reduction: rho = tanh(eta), g = cosh(eta)^{1+lam} chi(eta) gives
    [-Delta_{H^n} - (n-1)^2/4 + W] chi = -kappa^2 chi, kappa = lam-(d-1)/2,
    W = V(tanh eta)/cosh^2 eta   (checked numerically at random points).
 7. SUSY partner potentials: Vtilde = -4(d-1)/(d-2+rho^2) in dimension d+4,
    spectral parameter lam+1, d = 3..8.
 8. d = 3, lam = -2: explicit smooth eigenfunction (u-picture)
       phi(rho) = [3(1+rho^2)^3 arctan(rho) - 3 rho - 4 rho^3 - 9 rho^5] / (rho^2 (1+rho^2)).
"""
import sympy as sp

rho, eta, lam, t, T, r = sp.symbols("rho eta lambda t T r", positive=True)


def check_profiles():
    for d in range(3, 11):
        f0 = 2 * sp.atan(rho / sp.sqrt(d - 2))
        ode = ((rho**2 - 1) * sp.diff(f0, rho, 2) + (2 * rho - (d - 1) / rho) * sp.diff(f0, rho)
               + (d - 1) * sp.sin(2 * f0) / (2 * rho**2))
        print(f"1. profile d={d}:", sp.simplify(sp.expand_trig(ode)))
    # derivation of the self-similar ODE itself
    d = sp.symbols("d", positive=True)
    f = sp.Function("f")
    u = f(r / (T - t))
    pde = sp.diff(u, t, 2) - sp.diff(u, r, 2) - (d - 1) / r * sp.diff(u, r) + (d - 1) * sp.sin(2 * u) / (2 * r**2)
    pde = sp.simplify((pde * (T - t) ** 2).subs(r, rho * (T - t)))
    F = f(rho)
    ode = ((rho**2 - 1) * sp.diff(F, rho, 2) + (2 * rho - (d - 1) / rho) * sp.diff(F, rho)
           + (d - 1) * sp.sin(2 * F) / (2 * rho**2))
    print("1'. self-similar ODE derivation:", sp.simplify(pde - ode))


def check_reduction_and_modes():
    d = sp.symbols("d", positive=True)
    n = d + 2
    psi = sp.Function("psi")
    U = r * psi(t, r)
    expr = sp.diff(U, t, 2) - sp.diff(U, r, 2) - (d - 1) / r * sp.diff(U, r) + (d - 1) * sp.sin(2 * U) / (2 * r**2)
    target = r * (sp.diff(psi(t, r), t, 2) - sp.diff(psi(t, r), r, 2) - (n - 1) / r * sp.diff(psi(t, r), r)
                  + (d - 1) * (sp.sin(2 * r * psi(t, r)) - 2 * r * psi(t, r)) / (2 * r**3))
    print("2. u = r psi reduction:", sp.simplify(sp.expand(expr - target)))
    g, V = sp.Function("g"), sp.Function("V")
    w = (T - t) ** (-1 - lam) * g(r / (T - t))
    pde = sp.diff(w, t, 2) - sp.diff(w, r, 2) - (n - 1) / r * sp.diff(w, r) + V(r / (T - t)) / (T - t) ** 2 * w
    pde = sp.simplify((pde * (T - t) ** (3 + lam)).subs(r, rho * (T - t)))
    G = g(rho)
    ode = ((1 - rho**2) * sp.diff(G, rho, 2) + ((n - 1) / rho - 2 * (lam + 2) * rho) * sp.diff(G, rho)
           - (lam + 1) * (lam + 2) * G - V(rho) * G)
    print("3. mode equation:", sp.simplify(pde + ode))


def check_potential_gauge():
    for d in range(3, 9):
        n = d + 2
        f0 = 2 * sp.atan(rho / sp.sqrt(d - 2))
        V = (d - 1) * (sp.cos(2 * f0) - 1) / rho**2
        Vc = -8 * (d - 1) * (d - 2) / (d - 2 + rho**2) ** 2
        print(f"4. potential d={d}:", sp.simplify(sp.expand_trig(V) - Vc))
        g = 1 / (d - 2 + rho**2)
        L = 1
        res = ((1 - rho**2) * sp.diff(g, rho, 2) + ((n - 1) / rho - 2 * (L + 2) * rho) * sp.diff(g, rho)
               - (L + 1) * (L + 2) * g - Vc * g)
        print(f"5. gauge mode d={d}:", sp.simplify(res))


def check_hyperbolic():
    import random
    chi = sp.Function("chi")
    for (lv, d, e) in [(0.37, 3, 0.8), (1.3 + 0.4j, 4, 1.7), (-0.6, 5, 0.3), (2.2 - 1.1j, 7, 2.4)]:
        n = d + 2
        kap = lv - (d - 1) / 2
        Vf = lambda x: -8 * (d - 1) * (d - 2) / (d - 2 + x**2) ** 2
        test = sp.exp(-eta) * sp.cos(eta) + eta**2  # arbitrary test function
        gg = sp.cosh(eta) ** (1 + lv) * test
        D = lambda f: sp.cosh(eta) ** 2 * sp.diff(f, eta)  # d/drho
        rr = sp.tanh(eta)
        lhs = ((1 - rr**2) * D(D(gg)) + ((n - 1) / rr - 2 * (lv + 2) * rr) * D(gg)
               - (lv + 1) * (lv + 2) * gg - Vf(rr) * gg)
        hyp = (sp.diff(test, eta, 2) + (n - 1) / sp.tanh(eta) * sp.diff(test, eta)
               + ((n - 1) ** 2 / 4 - kap**2) * test - Vf(rr) / sp.cosh(eta) ** 2 * test)
        val = complex(sp.N((lhs / sp.cosh(eta) ** (lv + 3) - hyp).subs(eta, e)))
        print(f"6. hyperbolic reduction (lam={lv}, d={d}): |residual| = {abs(val):.1e}")


def check_susy():
    for d in range(3, 9):
        n = d + 2
        W = (-8 * (d - 1) * (d - 2) / (d - 2 + sp.tanh(eta) ** 2) ** 2) / sp.cosh(eta) ** 2
        m = sp.Rational(n - 1, 2)
        Q = m * (m - 1) / sp.sinh(eta) ** 2 + W
        kap1 = 1 - sp.Rational(d - 1, 2)
        psi1 = sp.sinh(eta) ** m / ((d - 2) * sp.cosh(eta) ** 2 + sp.sinh(eta) ** 2)
        r1 = sp.simplify((-sp.diff(psi1, eta, 2) + (Q + kap1**2) * psi1).rewrite(sp.exp))
        wsp = -sp.diff(psi1, eta) / psi1
        Qt = Q + 2 * sp.diff(wsp, eta)
        # partner written as a mode problem in dimension n+2: centrifugal (m+1)m/sinh^2
        Vt = sp.cosh(eta) ** 2 * (Qt - (m + 1) * m / sp.sinh(eta) ** 2)
        Vt_claim = -4 * (d - 1) / (d - 2 + sp.tanh(eta) ** 2)
        r2 = sp.simplify((Vt - Vt_claim).rewrite(sp.exp))
        print(f"7. SUSY d={d}: H psi1 = -kappa1^2 psi1 residual {r1}; partner potential residual {r2}")


def check_lambda_minus2():
    lamv = -2
    c2f = (1 - 6 * rho**2 + rho**4) / (1 + rho**2) ** 2  # cos(2 f0), d = 3
    phi = (3 * (1 + rho**2) ** 3 * sp.atan(rho) - 3 * rho - 4 * rho**3 - 9 * rho**5) / (rho**2 * (1 + rho**2))
    # u-picture (s = 0, R^3): v = (T-t)^{-lam} phi(rho)
    res = ((1 - rho**2) * sp.diff(phi, rho, 2) + (2 / rho - 2 * (lamv + 1) * rho) * sp.diff(phi, rho)
           - lamv * (lamv + 1) * phi - 2 * c2f / rho**2 * phi)
    print("8. d=3, lam=-2 eigenfunction residual:", sp.simplify(res))
    print("   Taylor at 0:", sp.series(phi, rho, 0, 6))
    print("   value/derivative at rho=1:", sp.N(phi.subs(rho, 1)), sp.N(sp.diff(phi, rho).subs(rho, 1)))
    # the cross-check that the same lam = -2 is also a zero of the Wronskian is in
    # geom_wavemaps_spectrum.py; here also check cos(2 f0) formula
    f0 = 2 * sp.atan(rho)
    print("   cos(2f0) formula:", sp.simplify(sp.expand_trig(sp.cos(2 * f0)) - c2f))


if __name__ == "__main__":
    check_profiles()
    check_reduction_and_modes()
    check_potential_gauge()
    check_hyperbolic()
    check_susy()
    check_lambda_minus2()
