#!/usr/bin/env python3
"""Symbolic checks for Section 5 (power family, Klein--Gordon, lower-order terms).

Run:  python3 power_symbolic.py
Every check prints 'OK' or raises AssertionError.

Checks
 1. Similarity-coordinate system for u_tt + gamma u_t + m^2 u - Delta u = F(u)
    (d = 2, general T):  extra terms -gamma*rho*psi2 - m^2 rho^2 psi1, rho = T e^{-tau}.
 2. ODE blowup: kappa and V0 = F'(u_T)(T-t)^2 = (s+1)(s+2) for
    F = |u|^{p-1}u (p>1) and for the MEMS force F = -lam w^{-q} (q>1).
 3. l-sector reduction of the eigenvalue PDE to Gauss' hypergeometric equation
    (d = 3, l = 0,1,2 with explicit solid harmonics; d = 2, l = 0..2).
 4. Polynomial eigenfunctions for alpha = -n / beta = -n solve the PDE.
 5. Constant (spatially homogeneous) sector: 2x2 matrix, eigenvalues 1, -2s-2.
 6. Jordan blocks at the collision exponents s = -1/2 (lambda=-1) and s = 1/2
    (lambda = -3), l = 0.
 7. Expansion of the spatially homogeneous blowup solution of
    u'' + gamma u' + m^2 u = F(u):  u = kappa rho^{-s}(1 + a1 rho + a2 rho^2 + ...),
    resonance at j = 2s+2, log terms (p = 3 with damping; MEMS q = 3 with damping).
 8. Explicit non-ODE profiles (cubic, all d >= 5; quadratic, d = 9).
 9. Tilted ("boosted") plane-wave family u = W(T - t + a.x) for the perturbed equation.
"""
import sympy as sp

ok = lambda msg: print(f"OK  {msg}")

# ---------------------------------------------------------------- 1
def check_similarity_system():
    """Identity check with an explicit (generic) psi1: for u = rho^{-s} psi1(tau, x/rho),
    psi2 := rho^{s+1} u_t satisfies
       d_tau psi1 = -Lambda psi1 - s psi1 + psi2,
       d_tau psi2 = Lap psi1 - Lambda psi2 - (s+1) psi2
                    + rho^{s+2} (u_tt - Lap u)          (exact identity),
    hence, if u_tt + gamma u_t + m^2 u - Lap u = F(u),
       d_tau psi2 = ... - gamma rho psi2 - m^2 rho^2 psi1 + rho^{s+2} F(u),
    with rho = T - t = T e^{-tau}."""
    t, x, y, T = sp.symbols('t x y T', positive=True)
    s = sp.Rational(3, 7)          # generic non-integer exponent
    tau, X, Y = sp.symbols('tau X Y', real=True)
    psi1 = sp.exp(sp.Rational(1, 3)*tau)*X**2*Y + sp.sin(tau + X)*Y**3 + sp.cos(X*Y)
    rho = T - t
    tau_t = -sp.log(T - t) + sp.log(T)
    phys = {tau: tau_t, X: x/rho, Y: y/rho}
    u = rho**(-s)*psi1.subs(phys)
    psi2_phys = rho**(s+1)*sp.diff(u, t)
    sim = {t: T - T*sp.exp(-tau), x: X*T*sp.exp(-tau), y: Y*T*sp.exp(-tau)}
    Lam = lambda g: X*sp.diff(g, X) + Y*sp.diff(g, Y)
    psi2 = sp.simplify(psi2_phys.subs(sim))
    # first equation
    e1 = sp.simplify(sp.diff(psi1, tau) - (-Lam(psi1) - s*psi1 + psi2))
    assert e1 == 0, e1
    # d_tau at fixed xi  ==  rho d_t - x.grad_x  in physical variables
    dtau_psi2 = (rho*sp.diff(psi2_phys, t) - x*sp.diff(psi2_phys, x) - y*sp.diff(psi2_phys, y)).subs(sim)
    box_u = (sp.diff(u, t, 2) - sp.diff(u, x, 2) - sp.diff(u, y, 2)).subs(sim)
    rhoS = T*sp.exp(-tau)
    rhs = (sp.diff(psi1, X, 2) + sp.diff(psi1, Y, 2) - Lam(psi2) - (s+1)*psi2 + rhoS**(s+2)*box_u)
    e2 = sp.simplify(dtau_psi2 - rhs)
    assert e2 == 0, e2
    # lower-order terms: rho^{s+2} gamma u_t = gamma rho psi2 ; rho^{s+2} m^2 u = m^2 rho^2 psi1
    ut = sp.diff(u, t).subs(sim); uu = u.subs(sim)
    assert sp.simplify(rhoS**(s+2)*ut - rhoS*psi2) == 0
    assert sp.simplify(rhoS**(s+2)*uu - rhoS**2*psi1) == 0
    ok("similarity system: d_tau psi2 = Lap psi1 - Lambda psi2 - (s+1)psi2 - gamma T e^{-tau} psi2"
       " - m^2 T^2 e^{-2tau} psi1 + rho^{s+2}F(u)")

# ---------------------------------------------------------------- 2
def check_ode_blowup():
    rho, p, q, lam = sp.symbols('rho p q lambda', positive=True)
    # power case
    s = 2/(p-1)
    kap = (s*(s+1))**(1/(p-1))
    uT = kap*rho**(-s)
    res = sp.simplify(sp.diff(uT, rho, 2) - uT**p)
    assert sp.simplify(res.subs(p, 3)) == 0 and sp.simplify(res.subs(p, sp.Rational(7, 3))) == 0
    V0 = sp.simplify(sp.powdenest(p*kap**(p-1), force=True))      # F'(u_T) rho^2
    assert sp.simplify(V0 - (s+1)*(s+2)) == 0
    ok("power case: kappa^(p-1) = s(s+1), V0 = p s(s+1) = (s+1)(s+2)")
    # MEMS case  w_tt - Lap w = -lam w^{-q}, w = kappa rho^sigma
    sig = 2/(q+1)
    kapM = (lam/(sig*(1-sig)))**(1/(q+1))
    w = kapM*rho**sig
    resM = sp.simplify(sp.diff(w, rho, 2) + lam*w**(-q))
    assert sp.simplify(resM.subs(q, 2)) == 0 and sp.simplify(resM.subs(q, 5)) == 0
    V0M = sp.simplify(sp.powdenest(sp.expand_power_base((lam*q*w**(-q-1))*rho**2, force=True), force=True))    # F'(w) = +lam q w^{-q-1}
    sM = -sig
    assert sp.simplify(V0M - (sM+1)*(sM+2)) == 0
    ok("MEMS case: s=-sigma, kappa^(q+1)=lam/(sigma(1-sigma)), V0 = q sigma(1-sigma) = (1-sigma)(2-sigma) = (s+1)(s+2)")
    # unified sign: F_p(u)=sgn(p-1)|u|^{p-1}u  => kappa^{p-1} = sgn(p-1) s(s+1)
    for pv in [sp.Rational(-3), sp.Rational(-2), sp.Rational(-3, 2), sp.Rational(1, 2), sp.Rational(-1, 2)]:
        sv = 2/(pv-1)
        print(f"    p={pv}: s={sv}, s(s+1)*sgn(p-1) = {sv*(sv+1)*sp.sign(pv-1)}  (needs >0 for ODE blowup profile)")

# ---------------------------------------------------------------- 3,4
def eig_pde_operator(g, coords, a, V0):
    """(delta_ij - xi_i xi_j) d_i d_j g - 2(a+1) Lambda g - (a(a+1)-V0) g"""
    d = len(coords)
    out = 0
    for i, xi in enumerate(coords):
        for j, xj in enumerate(coords):
            out += ((1 if i == j else 0) - xi*xj)*sp.diff(g, xi, xj)
    Lam = sum(xi*sp.diff(g, xi) for xi in coords)
    return out - 2*(a+1)*Lam - (a*(a+1)-V0)*g

def check_hypergeometric_reduction():
    a, V0 = sp.symbols('a V0')
    z = sp.symbols('z')
    h = sp.Function('h')
    x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
    harmonics3 = {0: sp.Integer(1), 1: x3, 2: x1*x2, 3: x1*x2*x3}
    harmonics2 = {0: sp.Integer(1), 1: x1, 2: x1*x2}
    for coords, harm in [((x1, x2, x3), harmonics3), ((x1, x2), harmonics2)]:
        d = len(coords)
        r2 = sum(c**2 for c in coords)
        for l, P in harm.items():
            g = P*h(r2)
            lhs = sp.expand(eig_pde_operator(g, coords, a, V0))
            A = a + l
            c = sp.Rational(d, 2) + l
            hyp = (z*(1-z)*sp.diff(h(z), z, 2) + (c - (A + sp.Rational(3, 2))*z)*sp.diff(h(z), z)
                   - sp.Rational(1, 4)*(A*(A+1) - V0)*h(z))
            rhs = 4*P*hyp.subs(z, r2)
            dd = sp.simplify((lhs - sp.expand(rhs.doit())).doit())
            assert dd == 0, (d, l, dd)
    ok("l-sector reduction: g = r^l Y_l h(r^2) -> z(1-z)h'' + [c-(A+3/2)z]h' - (A(A+1)-V0)/4 h = 0, "
       "c = d/2+l, A = a+l   (checked d=3: l=0..3, d=2: l=0..2)")

    # eigenfunctions
    s_ = sp.Rational(3, 7)
    V0v = (s_+1)*(s_+2)
    for coords, harm in [((x1, x2, x3), harmonics3), ((x1, x2), harmonics2)]:
        d = len(coords)
        r2 = sum(cc**2 for cc in coords)
        for l, P in harm.items():
            for n in range(3):
                for fam in (1, 2):
                    lamv = (1 - l - 2*n) if fam == 1 else (-2*s_ - 2 - l - 2*n)
                    av = s_ + lamv
                    A = av + l
                    al = (A - (s_+1))/2
                    be = (A + (s_+2))/2
                    if fam == 1:
                        assert al == -n
                    else:
                        assert be == -n
                    c = sp.Rational(d, 2) + l
                    hpoly = sp.hyper([al, be], [c], z)
                    hpoly = sp.expand(sp.hyperexpand(hpoly))
                    g = P*hpoly.subs(z, r2)
                    resid = sp.simplify(sp.expand(eig_pde_operator(g, coords, av, V0v)))
                    assert resid == 0, (d, l, n, fam, resid)
    ok("eigenfunctions r^l Y_l 2F1(alpha,beta;d/2+l;r^2) with alpha=-n (lambda=1-l-2n) "
       "or beta=-n (lambda=-2s-2-l-2n) solve the eigen-PDE (s=3/7, d=2,3)")

# ---------------------------------------------------------------- 5
def check_constant_sector():
    s, lam = sp.symbols('s lambda')
    V0 = (s+1)*(s+2)
    M = sp.Matrix([[-s, 1], [V0, -(s+1)]])
    ev = M.eigenvals()
    assert set(sp.simplify(e) for e in ev) == {1, -2*s-2}, ev
    ok("constant sector: L restricted to constants = [[-s,1],[V0,-(s+1)]], eigenvalues {1, -2s-2}")

# ---------------------------------------------------------------- 6
def L_apply(phi1, phi2, coords, s, V0):
    Lam = lambda g: sum(c*sp.diff(g, c) for c in coords)
    Lap = lambda g: sum(sp.diff(g, c, 2) for c in coords)
    return (-Lam(phi1) - s*phi1 + phi2,
            Lap(phi1) - Lam(phi2) - (s+1)*phi2 + V0*phi1)

def check_jordan(s_val, lam_val, d=3, deg=6, l=0):
    """Is there a generalized eigenvector (L - lam) phi = e for the eigenvector e
    in the l-sector (l = 0: radial, l = 1: prefactor x1)?  Returns (g, solutions)."""
    coords = sp.symbols('x1:%d' % (d+1), real=True)
    r2 = sum(c**2 for c in coords)
    P = sp.Integer(1) if l == 0 else coords[0]
    V0 = (s_val+1)*(s_val+2)
    z = sp.symbols('z')
    A = s_val + lam_val + l
    al, be = (A - (s_val+1))/2, (A + (s_val+2))/2
    hpoly = sp.expand(sp.hyperexpand(sp.hyper([al, be], [sp.Rational(d, 2) + l], z)))
    g = P*hpoly.subs(z, r2)
    Lam = lambda f: sum(c*sp.diff(f, c) for c in coords)
    e1, e2 = g, sp.expand((lam_val + s_val)*g + Lam(g))
    r1, r2_ = L_apply(e1, e2, coords, s_val, V0)
    assert sp.simplify(r1 - lam_val*e1) == 0 and sp.simplify(r2_ - lam_val*e2) == 0
    cs = sp.symbols('c0:%d' % (deg+1))
    phi1 = P*sum(cs[k]*r2**k for k in range(deg+1))
    phi2 = e1 + Lam(phi1) + (s_val + lam_val)*phi1
    o1, o2 = L_apply(phi1, phi2, coords, s_val, V0)
    eq = sp.expand(o2 - lam_val*phi2 - e2)
    rr = sp.symbols('rr')
    eqr = sp.expand(eq.subs({coords[0]: rr, **{c: 0 for c in coords[1:]}}))
    sol = sp.solve(sp.Poly(eqr, rr).coeffs(), cs, dict=True)
    return g, sol


def check_jordan_blocks():
    for s_val, lam_val in [(sp.Rational(-1, 2), -1), (sp.Rational(1, 2), -3)]:
        g, sol = check_jordan(s_val, lam_val)
        assert sol, "no generalized eigenvector"
        print(f"    s={s_val}, lambda={lam_val}: eigenfunction g={g}, generalized eigvec phi1 coeffs {sol[0]}")
    # controls: no Jordan block at non-collision eigenvalues
    for s_val in (sp.Rational(1, 3), sp.Rational(-2, 3), sp.Integer(1), sp.Rational(-1, 2)):
        assert not check_jordan(s_val, 1)[1]            # time translation, l = 0
        assert not check_jordan(s_val, 0, l=1)[1]       # Lorentz boosts, l = 1
    assert not check_jordan(sp.Rational(1, 3), -1)[1]
    assert not check_jordan(sp.Rational(-2, 3), sp.Rational(-2, 3))[1]   # MEMS q=2 gap eigenvalue
    ok("Jordan block (algebraic mult. 2, geometric 1) at the collision exponents "
       "s=-1/2 (lambda=-1, MEMS q=3) and s=1/2 (lambda=-3, p=5); lambda=1 (l=0) and "
       "lambda=0 (l=1) semisimple for s in {1/3,-2/3,1,-1/2}; lambda=-2/3 simple for q=2")


# ---------------------------------------------------------------- 7
def ode_expansion(pv, gam, m, J=5, sign=1, extra_log=None):
    """u'' (in rho) - gam u_rho + m^2 u = sign*u^p,  u = kappa rho^{-s} sum a_j rho^j.
    returns the list of a_j (sympy), stops at resonance j=2s+2 if integer."""
    s = 2/(pv-1)
    rho = sp.symbols('rho', positive=True)
    kap = (sign*s*(s+1))**(1/(pv-1))
    a = sp.symbols('a1:%d' % (J+1))
    series = 1 + sum(a[j-1]*rho**j for j in range(1, J+1))
    u = kap*rho**(-s)*series
    lhs = sp.diff(u, rho, 2) - gam*sp.diff(u, rho) + m**2*u
    # sign*u^p = sign*kap^p rho^{-s p} series^p ; s p = s + 2
    rhs = sign*kap**pv*rho**(-s-2)*sp.series(series**pv, rho, 0, J+1).removeO()
    expr = sp.expand(sp.simplify((lhs - rhs)*rho**(s+2)/kap))
    sols = {}
    for j in range(1, J+1):
        cj = sp.expand(expr.coeff(rho, j).subs(sols))
        if sp.simplify(sp.diff(cj, a[j-1])) == 0:
            sols['resonance'] = (j, sp.simplify(cj))
            break
        sols[a[j-1]] = sp.solve(cj, a[j-1])[0]
    return s, sols

def check_ode_expansion():
    gam, m = sp.symbols('gamma m', positive=True)
    pp = sp.symbols('p', positive=True)
    # general p: a1, a2 in closed form via generic s
    s = sp.symbols('s', positive=True)
    rho = sp.symbols('rho', positive=True)
    a1, a2 = sp.symbols('a1 a2')
    # generic order-by-order recursion: coefficient of a_j is (j-s)(j-s-1)-(s+1)(s+2)
    j = sp.symbols('j')
    D = sp.factor((j-s)*(j-s-1) - (s+1)*(s+2))
    print("    order-j operator:", D, " -> resonances j=-1 (time translation), j=2s+2")
    for pv in [3, 5, 7, sp.Rational(7, 3)]:
        sv, sols = ode_expansion(sp.Integer(pv) if not isinstance(pv, sp.Rational) else pv, gam, m, J=4)
        print(f"    p={pv}, s={sv}:", {str(k): sp.factor(v) if not isinstance(v, tuple) else v for k, v in sols.items()})
    # general formula check:
    for pv in [3, 5, 7, sp.Rational(7, 3), sp.Rational(9, 2)]:
        sv, sols = ode_expansion(sp.nsimplify(pv), gam, m, J=2)
        A1 = gam*sv/(2*(2*sv+1))
        assert sp.simplify(sols[sp.Symbol('a1')] - A1) == 0
    ok("a1 = gamma s / (2(2s+1))  (power family, all tested p)")
    # MEMS q=2 and q=3 (F = -w^{-q}, i.e. p = -q with sign = -1)
    for qv in [2, 3, 5]:
        sv, sols = ode_expansion(sp.Integer(-qv), gam, m, J=3, sign=-1)
        print(f"    MEMS q={qv}, s={sv}:", {str(k): sp.factor(v) if not isinstance(v, tuple) else v for k, v in sols.items()})


def check_a2_general():
    """closed forms for a1, a2 with symbolic s (u^p expanded with p = 1+2/s)."""
    s, gam, m = sp.symbols('s gamma m', positive=True)
    p = 1 + 2/s
    a1, a2 = sp.symbols('a1 a2')
    # equation for coefficients: sum over j of [(j-s)(j-s-1) a_j + gam (s-j+1) a_{j-1} + m^2 a_{j-2}]
    #   = s(s+1) [ (1+a1 r + a2 r^2)^p ]_j
    # binomial: (1+X)^p = 1 + pX + p(p-1)/2 X^2
    c1 = (1-s)*(-s)*a1 + gam*s - s*(s+1)*p*a1
    sol1 = sp.solve(c1, a1)[0]
    c2 = (2-s)*(1-s)*a2 + gam*(s-1)*sol1 + m**2 - s*(s+1)*(p*a2 + p*(p-1)/2*sol1**2)
    sol2 = sp.factor(sp.solve(c2, a2)[0])
    print("    a1 =", sp.factor(sol1))
    print("    a2 =", sol2)
    # numeric cross-check with ode_expansion for p=3
    _, sols = ode_expansion(sp.Integer(3), gam, m, J=2)
    assert sp.simplify(sols[sp.Symbol('a2')] - sol2.subs(s, 1)) == 0
    _, sols = ode_expansion(sp.Integer(7), gam, m, J=2)
    assert sp.simplify(sols[sp.Symbol('a2')] - sol2.subs(s, sp.Rational(1, 3))) == 0
    ok("closed forms of a1, a2 (checked against direct expansion for p=3,7)")
    return sol1, sol2

# ---------------------------------------------------------------- 8
def check_explicit_profiles():
    """Non-ODE self-similar profiles: cubic (Glogic-Schoerkhuber form) and quadratic d=9."""
    r = sp.symbols('r', positive=True)
    d = sp.symbols('d', positive=True)
    def resid(U, s, p, dd):
        Lam = lambda f: r*sp.diff(f, r)
        Lap = lambda f: sp.diff(f, r, 2) + (dd-1)/r*sp.diff(f, r)
        A = s*U + Lam(U)
        return sp.simplify((s+1)*A + Lam(A) - Lap(U) - U**p)
    Uc = 2*sp.sqrt(2)*sp.sqrt((d-1)*(d-4))/(d - 4 + 3*r**2)
    assert resid(Uc, 1, 3, d) == 0
    Uq = 240*(7 - r**2)/(7 + 5*r**2)**2
    assert resid(Uq, 2, 2, 9) == 0
    ok("explicit profiles: u=(T-t)^{-1} 2sqrt2 sqrt((d-1)(d-4))/(d-4+3rho^2) solves u_tt-Lap u=u^3 (all d); "
       "u=(T-t)^{-2} 240(7-rho^2)/(7+5rho^2)^2 solves u_tt-Lap u=u^2 in d=9")

def check_tilted_family():
    """u = W(T - t + a.x) reduces u_tt + gam u_t + m^2 u - Lap u = F(u) to the ODE
    (1-|a|^2) W'' - gam W' + m^2 W = F(W); for gam = m = 0, W = kappa (1-|a|^2)^{s/2} phi^{-s}."""
    t, x, y, T, a1, a2, gam, m = sp.symbols('t x y T a1 a2 gamma m', real=True)
    W = sp.Function('W')
    phi = T - t + a1*x + a2*y
    u = W(phi)
    expr = sp.diff(u, t, 2) + gam*sp.diff(u, t) + m**2*u - sp.diff(u, x, 2) - sp.diff(u, y, 2)
    ph = sp.symbols('phi')
    target = (1 - a1**2 - a2**2)*sp.Derivative(W(ph), (ph, 2)) - gam*sp.Derivative(W(ph), ph) + m**2*W(ph)
    # compare via substitution of a concrete W
    Wc = lambda z: sp.exp(sp.Rational(1, 3)*z) + z**5
    uc = Wc(phi)
    e1 = sp.diff(uc, t, 2) + gam*sp.diff(uc, t) + m**2*uc - sp.diff(uc, x, 2) - sp.diff(uc, y, 2)
    zz = sp.symbols('zz')
    e2 = ((1 - a1**2 - a2**2)*sp.diff(Wc(zz), zz, 2) - gam*sp.diff(Wc(zz), zz) + m**2*Wc(zz)).subs(zz, phi)
    assert sp.simplify(e1 - e2) == 0
    p = sp.Integer(5); s = 2/(p-1); kap = (s*(s+1))**(1/(p-1))
    A2 = a1**2 + a2**2
    Wb = kap*(1 - A2)**(s/2)*ph**(-s)
    assert sp.simplify((1 - A2)*sp.diff(Wb, ph, 2) - Wb**p) == 0
    ok("tilted family u = W(T-t+a.x): ODE (1-|a|^2)W'' - gamma W' + m^2 W = F(W); "
       "exact boosted ODE profile kappa(1-|a|^2)^{s/2}(T-t+a.x)^{-s} (checked p=5)")

def check_q3_log():
    """MEMS q=3 (s=-1/2) with damping: U = sqrt2 rho^{1/2}(1 + b rho log rho + c rho + ...), b = gamma/4."""
    rho, gam, m, b, c, L = sp.symbols('rho gamma m b c L', positive=True)
    U = sp.sqrt(2)*rho**sp.Rational(1, 2)*(1 + b*rho*sp.log(rho) + c*rho)
    expr = sp.diff(U, rho, 2) - gam*sp.diff(U, rho) + m**2*U + U**(-3)
    e = sp.expand(sp.series(sp.simplify(expr*rho**sp.Rational(3, 2)/sp.sqrt(2)), rho, 0, 2).removeO())
    e = sp.expand(e.subs(sp.log(rho), L))
    assert sp.solve(e.coeff(rho, 1).coeff(L, 0), b) == [gam/4]
    ok("MEMS q=3 damped ODE blowup: log coefficient b = gamma/4 (c free)")


if __name__ == "__main__":
    check_similarity_system()
    check_ode_blowup()
    check_hypergeometric_reduction()
    check_constant_sector()
    check_jordan_blocks()
    check_ode_expansion()
    check_a2_general()
    check_explicit_profiles()
    check_tilted_family()
    check_q3_log()
