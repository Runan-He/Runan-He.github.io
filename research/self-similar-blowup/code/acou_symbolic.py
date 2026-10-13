#!/usr/bin/env python3
"""Symbolic checks for Section 7 (nonlinear acoustics, derivative nonlinearities).

Run:  python3 acou_symbolic.py      (prints 'OK ...' for every check, raises otherwise)

Model (b):  u_tt - Lap u = (u_t)^2        ("Riccati" / semilinear Westervelt)
Model (c):  u_tt - Lap u = (u_t)^2 - |grad u|^2   (Nirenberg null form)
Westervelt: (1-2k p) p_tt - c^2 Lap p - b Lap p_t = 2k (p_t)^2
Kuznetsov:  psi_tt - c^2 Lap psi - b Lap psi_t = d_t( K psi_t^2 + |grad psi|^2 ),  K = (B/2A)/c^2

Checks
 B1  ODE blowup u = -log(T-t) + c; tilted family u_a = -(1-|a|^2) log(T-t+a.x) (d=3);
     similarity-variable form of u_a: psi1 = -|a|^2 tau + ..., psi2 = (1-|a|^2)/(1+a.xi).
 B2  similarity system (d=2, generic psi1): d_tau psi1 = -Lam psi1 + psi2 - 1,
     d_tau psi2 = Lap psi1 - Lam psi2 - psi2 + psi2^2 + rho^2 (u_tt - Lap u - u_t^2).
 B3  linearised operator L(phi) = (-Lam phi1 + phi2, Lap phi1 - Lam phi2 + phi2) is
     conjugate to the free operator L0^{(s=-1)} via S(P1,P2) = (P1, P1+P2);
     physical version: v_tt - Lap v - (2/rho) v_t = 0  <=>  w = rho v is a free wave.
 B4  eigenvalue ODE -> hypergeometric with alpha=(l+lam)/2, beta=(l+lam-1)/2, c=l+d/2;
     polynomial eigenfunctions from homogeneous polynomial free waves of degree j,
     lambda = 1 - j, for d = 1, 3 and j <= 3.
 B5  no Jordan chains for lambda = 1 (l=0) and lambda = 0 (l=0,1): the equation
     (L - lam) f = g_lam has only solutions with log / non-smooth behaviour at r = 1.
 B6  quadratic drift: the projection mechanism; check that psi(u_a) solves the
     nonlinear system and that its O(a^2) part contains -|a|^2 tau (shift direction).
 B7  viscous term b Lap u_t in similarity variables equals b e^{tau} Lap psi2 / T (dominant).
 C1  Nirenberg transform: Box(e^{-u}) = e^{-u}(u_t^2 - |grad u|^2 - Box u);
     null-form tilted family -log(T-t+a.x) has coefficient 1 (Lorentz invariant);
     -log((T-t)^2 + |x|^2/d) is an exact (non-ODE) self-similar blowup solution.
 C2  model (b) in the variable w = e^{-u}:  w Box w = -|grad_x w|^2.
 E1  energy identity for the free system with exponent s' (d = 2, pointwise):
     2 grad f1 . grad(-Lam f1 - s' f1 + f2) + 2 f2 (Lap f1 - Lam f2 - (s'+1) f2)
       = (d-2-2s')(|grad f1|^2 + f2^2) + div( -xi(|grad f1|^2+f2^2) + 2 f2 grad f1 ),
     and the commutator  D^beta Lam = (Lam + |beta|) D^beta.
 W1  Westervelt <=> (q^2)_tt = 2 c^2 Lap q for q = 1 - 2kp (inviscid); homogeneous
     solutions q^2 = A + B t; p_t ~ (T-t)^{-1/2}.
 W2  Kuznetsov - Westervelt (potential forms) = d_t(|grad psi|^2 - psi_t^2/c^2)
     (time derivative of the null form); Kuznetsov "uniform compression" solution
     psi = f(t) + B|x|^2/(2 d c^2) with (f' - K f'^2)' = B.
 W2b Kuznetsov in r = 1-2K psi_t: (r^2)_tt = 2c^2 Lap r - 8K d_t(grad psi . grad psi_t).
 W3  Westervelt self-similar variables eta = x/rho^{3/4}, q = rho^{1/2} Q:
     (d_tau + 3/4 D)(d_tau + 3/4 D - 1) Q^2 = 2 c^2 Lap_eta Q   (exact).
 W4  linearisation at Q0: modes e^{lam tau} h(zeta) with
     (Lam + 4lam/3)(Lam + 4(lam-1)/3) h = Lap h, hypergeometric parameters, explicit
     modes lam = 1, 1/4, 0, -1/2, -3/4.
"""
import sympy as sp

ok = lambda msg: print(f"OK  {msg}")


def box(u, t, xs):
    return sp.diff(u, t, 2) - sum(sp.diff(u, x, 2) for x in xs)


def grad2(u, xs):
    return sum(sp.diff(u, x)**2 for x in xs)


# ------------------------------------------------------------------ B1
def check_B1():
    t, T = sp.symbols('t T', real=True)
    xs = sp.symbols('x1:4', real=True)
    a = sp.symbols('a1:4', real=True)
    c0 = sp.Symbol('c0', real=True)
    uT = -sp.log(T - t) + c0
    assert sp.simplify(box(uT, t, xs) - sp.diff(uT, t)**2) == 0
    a2 = sum(ai**2 for ai in a)
    f = T - t + sum(ai*xi for ai, xi in zip(a, xs))
    ua = -(1 - a2)*sp.log(f) + c0
    res = sp.simplify(box(ua, t, xs) - sp.diff(ua, t)**2)
    assert res == 0, res
    # general coefficient: -C log f solves iff C = 1-|a|^2 (or C=0)
    C = sp.Symbol('C')
    uC = -C*sp.log(f)
    eq = sp.factor(sp.simplify((box(uC, t, xs) - sp.diff(uC, t)**2)*f**2))
    sol = sp.solve(eq, C)
    assert set(sol) == {0, 1 - a2} or set(sp.simplify(s_) for s_ in sol) == {0, sp.expand(1 - a2)}, sol
    # similarity variables, T = 1, rho = e^{-tau}
    tau = sp.Symbol('tau', real=True)
    xi = sp.symbols('xi1:4', real=True)
    rho = sp.exp(-tau)
    sub = {t: 1 - rho, T: 1, **{x: rho*q for x, q in zip(xs, xi)}}
    psi1 = sp.expand_log(sp.simplify((ua + sp.log(T - t)).subs(sub)), force=True)
    adot = sum(ai*q for ai, q in zip(a, xi))
    target1 = -a2*tau - (1 - a2)*sp.log(1 + adot) + c0
    assert sp.simplify(sp.expand_log(psi1 - target1, force=True)) == 0
    psi2 = sp.simplify(((T - t)*sp.diff(ua, t)).subs(sub))
    assert sp.simplify(psi2 - (1 - a2)/(1 + adot)) == 0
    ok("B1 ODE blowup -log(T-t)+c; tilted family -(1-|a|^2)log(T-t+a.x) (d=3) exact; "
       "psi1 = -|a|^2 tau - (1-|a|^2)log(1+a.xi) + c, psi2 = (1-|a|^2)/(1+a.xi)")


# ------------------------------------------------------------------ B2
def check_B2():
    """Exact identity with a generic psi1 (d = 2):
    u = psi1(tau, x/rho) - log rho, rho = T - t, tau = -log rho + log T."""
    t, x, y, T = sp.symbols('t x y T', positive=True)
    tau, X, Y = sp.symbols('tau X Y', real=True)
    psi1 = sp.exp(sp.Rational(1, 3)*tau)*X**2*Y + sp.sin(tau + X)*Y**3 + sp.cos(X*Y)
    rho = T - t
    phys = {tau: -sp.log(rho) + sp.log(T), X: x/rho, Y: y/rho}
    u = psi1.subs(phys) - sp.log(rho)
    sim = {t: T - T*sp.exp(-tau), x: X*T*sp.exp(-tau), y: Y*T*sp.exp(-tau)}
    Lam = lambda g: X*sp.diff(g, X) + Y*sp.diff(g, Y)
    psi2_phys = rho*sp.diff(u, t)
    psi2 = sp.simplify(psi2_phys.subs(sim))
    e1 = sp.simplify(sp.diff(psi1, tau) - (-Lam(psi1) + psi2 - 1))
    assert e1 == 0, e1
    dtau_psi2 = (rho*sp.diff(psi2_phys, t) - x*sp.diff(psi2_phys, x)
                 - y*sp.diff(psi2_phys, y)).subs(sim)
    R = (box(u, t, (x, y)) - sp.diff(u, t)**2).subs(sim)
    rhoS = T*sp.exp(-tau)
    rhs = (sp.diff(psi1, X, 2) + sp.diff(psi1, Y, 2) - Lam(psi2) - psi2 + psi2**2 + rhoS**2*R)
    e2 = sp.simplify(dtau_psi2 - rhs)
    assert e2 == 0, e2
    ok("B2 similarity system: d_tau psi1 = -Lam psi1 + psi2 - 1, "
       "d_tau psi2 = Lap psi1 - Lam psi2 - psi2 + psi2^2  (+ rho^2 * residual)")


# ------------------------------------------------------------------ B3
def check_B3():
    X, Y = sp.symbols('X Y', real=True)
    P1 = sp.Function('P1')(X, Y)
    P2 = sp.Function('P2')(X, Y)
    Lam = lambda g: X*sp.diff(g, X) + Y*sp.diff(g, Y)
    Lap = lambda g: sp.diff(g, X, 2) + sp.diff(g, Y, 2)
    L = lambda f1, f2: (-Lam(f1) + f2, Lap(f1) - Lam(f2) + f2)
    s = -1
    L0 = lambda f1, f2: (-Lam(f1) - s*f1 + f2, Lap(f1) - Lam(f2) - (s + 1)*f2)
    S = lambda f1, f2: (f1, f1 + f2)
    lhs = L(*S(P1, P2))
    rhs = S(*L0(P1, P2))
    assert all(sp.simplify(l - r) == 0 for l, r in zip(lhs, rhs))
    # physical version (d = 2)
    t, x, y, T = sp.symbols('t x y T', real=True)
    w = sp.Function('w')(t, x, y)
    rho = T - t
    v = w/rho
    expr = sp.simplify((box(v, t, (x, y)) - 2/rho*sp.diff(v, t))*rho - box(w, t, (x, y)))
    assert expr == 0, expr
    # linearisation of (u_t)^2 at u_T = -log(T-t): 2 u_T' v_t = 2 v_t/(T-t)
    uT = -sp.log(T - t)
    assert sp.simplify(2*sp.diff(uT, t) - 2/rho) == 0
    ok("B3 L = S L0^{(-1)} S^{-1};  v_tt - Lap v - (2/rho) v_t = rho^{-1} Box(rho v)")


# ------------------------------------------------------------------ B4
def check_B4():
    r, z, lam = sp.symbols('r z lambda')
    for d in (1, 2, 3, 5):
        for l in range(0, 4):
            if d == 1 and l > 1:
                continue
            h = sp.Function('h')
            g = r**l*h(r**2)
            # radial form of Lap - xi xi : D^2 - 2 lam Lam - lam(lam-1) on r^l Y_l h(r^2)
            gr = r**l*sp.Function('H')(r)
            Hf = sp.Function('H')(r)
            ode = ((1 - r**2)*sp.diff(gr, r, 2) + ((d - 1)/r - 2*lam*r)*sp.diff(gr, r)
                   - (l*(l + d - 2)/r**2 + lam*(lam - 1))*gr)
            ode = sp.simplify(ode/r**l)
            target = ((1 - r**2)*sp.diff(Hf, r, 2) + ((2*l + d - 1)/r - 2*(l + lam)*r)*sp.diff(Hf, r)
                      - (l + lam)*(l + lam - 1)*Hf)
            assert sp.simplify(ode - target) == 0
            # z = r^2 : hypergeometric with alpha=(l+lam)/2, beta=(l+lam-1)/2, c=l+d/2
            Gz = sp.Function('G')
            Hz = Gz(r**2)
            e = target.subs(Hf, Hz).doit()
            e = sp.simplify(e.subs(r, sp.sqrt(z)))
            al, be, cc = (l + lam)/2, (l + lam - 1)/2, l + sp.Rational(d, 2)
            G = Gz(z)
            hyp = z*(1 - z)*sp.diff(G, z, 2) + (cc - (al + be + 1)*z)*sp.diff(G, z) - al*be*G
            assert sp.simplify(e - 4*hyp) == 0, (d, l)
    ok("B4a eigen-ODE <-> 2F1(alpha,beta;c;r^2), alpha=(l+lam)/2, beta=(l+lam-1)/2, c=l+d/2")
    # polynomial modes: w homogeneous free wave of degree j  ->  lambda = 1 - j
    for d in (1, 3):
        xs = sp.symbols(f'x1:{d+1}', real=True)
        tt, T = sp.symbols('t T', real=True)
        rh = sp.Symbol('rho', positive=True)
        cands = {0: [sp.Integer(1)], 1: [rh] + list(xs),
                 2: [rh**2 + sum(x**2 for x in xs)/d, rh*xs[0]] + ([xs[0]*xs[1]] if d > 1 else []),
                 3: [rh**3 + 3*rh*sum(x**2 for x in xs)/d, xs[0]*(rh**2 + sum(x**2 for x in xs)/(d + 2))]}
        for j, ws in cands.items():
            for wp in ws:
                wt = wp.subs(rh, T - tt)
                assert sp.simplify(box(wt, tt, xs)) == 0, (d, j, wp)
                # eigenfunction in similarity variables: v = w/rho = rho^{j-1} g(xi)
                xi = sp.symbols(f'xi1:{d+1}', real=True)
                g1 = sp.simplify((wp/rh).subs({x: rh*q for x, q in zip(xs, xi)})/rh**(j - 1))
                assert not g1.has(rh), (d, j, g1)
                lamj = 1 - j
                Lam = lambda f: sum(q*sp.diff(f, q) for q in xi)
                Lap = lambda f: sum(sp.diff(f, q, 2) for q in xi)
                g2 = lamj*g1 + Lam(g1)
                assert sp.simplify(lamj*g2 - (Lap(g1) - Lam(g2) + g2)) == 0
    ok("B4b homogeneous polynomial free waves of degree j give eigenfunctions with lambda = 1-j "
       "(d = 1, 3; j <= 3; l = 0..j)")


# ------------------------------------------------------------------ B5
def check_B5():
    """Jordan chains. (L - lam) f = g with g the eigenfunction.  Elimination gives
    Lap f1 - xi xi:D^2 f1 - 2 lam Lam f1 - lam(lam-1) f1 = g2 + (Lam + lam - 1) g1.
    We solve the radial ODE exactly and show the solution regular at r = 0 is not
    smooth at r = 1 (log term), hence not in H^k for k >= k_d."""
    r = sp.Symbol('r', positive=True)
    out = []
    for d in (1, 3):
        for (lam, l, g1, g2) in [(1, 0, 1, 1), (0, 0, 1, 0), (0, 1, r, r)]:
            # f1 = r^l F(r)  (times Y_l); RHS also r^l * (...)
            F = sp.Function('F')
            Fr = F(r)
            lhs = ((1 - r**2)*sp.diff(Fr, r, 2) + ((2*l + d - 1)/r - 2*(l + lam)*r)*sp.diff(Fr, r)
                   - (l + lam)*(l + lam - 1)*Fr)
            # RHS: g2 + (Lam + lam - 1) g1, divided by r^l
            rhs = sp.simplify((g2 + r*sp.diff(g1, r) + (lam - 1)*g1)/r**l)
            # first-order ODE for F' when the zeroth-order coefficient vanishes
            coef0 = (l + lam)*(l + lam - 1)
            if coef0 == 0:
                G = sp.Function('G')
                odeG = sp.Eq((1 - r**2)*sp.diff(G(r), r) + ((2*l + d - 1)/r - 2*(l + lam)*r)*G(r), rhs)
                solG = sp.dsolve(odeG, G(r))
                Gs = solG.rhs
                # choose constant so that G is regular at r = 0 (and odd if d = 1, l = 0)
                C1 = sp.Symbol('C1')
                if (2*l + d - 1) > 0:
                    cond = sp.limit(sp.simplify(Gs*r**(2*l + d - 1)), r, 0)
                else:
                    cond = sp.limit(Gs, r, 0)       # parity: F' odd => F'(0) = 0
                cands = sp.solve([cond], C1, dict=True)
                Greg = sp.simplify(Gs.subs(cands[0]) if cands else Gs)
                assert not Greg.has(C1), Greg
                # behaviour near r = 1: look for log(eps) up to order eps^2
                eps = sp.Symbol('epsilon', positive=True)
                near1 = sp.expand(sp.series(Greg.subs(r, 1 - eps), eps, 0, 3).removeO())
                out.append((d, lam, l, sp.simplify(Greg), near1))
                neg_pow = any(isinstance(a_, sp.Pow) and a_.base == eps and a_.exp < 0
                              for a_ in sp.preorder_traversal(near1))
                assert near1.has(sp.log(eps)) or neg_pow, near1
            else:
                raise RuntimeError("unexpected")
    for d, lam, l, G, n1 in out:
        print(f"    d={d} lam={lam} l={l}:  F'(r) = {G}\n        near r=1-eps: {n1}")
    ok("B5 no Jordan chains for lam=1 (l=0), lam=0 (l=0,1) in d=1,3: generalized "
       "eigenfunction would need log(1-r) terms (not in H^k, k >= (d+3)/2)")


# ------------------------------------------------------------------ B6
def check_B6():
    """Second-order drift.  In similarity variables (T=1) phi = psi - (c0, 1).
    For u_a (d = 1, a small) expand psi1 in a: the a^2 term contains -a^2 tau
    (shift eigenvector (1,0), lambda = 0)."""
    tau, xi, a = sp.symbols('tau xi a', real=True)
    psi1 = -a**2*tau - (1 - a**2)*sp.log(1 + a*xi)
    psi2 = (1 - a**2)/(1 + a*xi)
    # exact solution of the nonlinear system:
    e1 = sp.simplify(sp.diff(psi1, tau) - (-xi*sp.diff(psi1, xi) + psi2 - 1))
    e2 = sp.simplify(sp.diff(psi2, tau) - (sp.diff(psi1, xi, 2) - xi*sp.diff(psi2, xi) - psi2 + psi2**2))
    assert e1 == 0 and e2 == 0
    s1 = sp.series(psi1, a, 0, 3).removeO()
    assert sp.simplify(s1.coeff(a, 1) + xi) == 0       # tilt mode -(xi, xi), lambda = 0
    assert sp.simplify(s1.coeff(a, 2) - (-tau + xi**2/2)) == 0
    # second-order equation for the a^2 coefficient: forced by N(phi^(1)) = (0, xi^2):
    # (L) applied to (q1, q2) = (-tau + xi^2/2, ...) ; check d_tau q1 = -1 comes from
    # projecting (0, xi^2) onto the shift direction.
    q1 = -tau + xi**2/2
    s2 = sp.series(psi2, a, 0, 3).removeO()
    q2 = s2.coeff(a, 2)
    r1 = sp.simplify(sp.diff(q1, tau) - (-xi*sp.diff(q1, xi) + q2))
    r2 = sp.simplify(sp.diff(q2, tau) - (sp.diff(q1, xi, 2) - xi*sp.diff(q2, xi) + q2 + xi**2))
    assert r1 == 0 and r2 == 0
    ok("B6 u_a solves the similarity system; O(a): eigenmode -(xi, xi) at lam=0; "
       "O(a^2): q1 = -tau + xi^2/2 solves d_tau q = L q + (0, xi^2): secular drift "
       "along the shift mode (1,0), not a Jordan block of L")


# ------------------------------------------------------------------ B7
def check_B7():
    t, x, T, b = sp.symbols('t x T b', positive=True)
    tau, X = sp.symbols('tau X', real=True)
    psi1 = sp.Function('p1')
    rho = T - t
    u = psi1(-sp.log(rho) + sp.log(T), x/rho) - sp.log(rho)
    term = b*sp.diff(u, x, 2, t)            # b Lap u_t (d = 1)
    sim = {t: T - T*sp.exp(-tau), x: X*T*sp.exp(-tau)}
    psi2 = (rho*sp.diff(u, t))
    lhs = sp.simplify((rho**2*term).subs(sim))
    rhs = sp.simplify((b/rho*sp.diff(psi2, x, 2)*rho**2).subs(sim))
    assert sp.simplify(lhs - rhs) == 0
    # i.e. rho^2 * b Lap_x u_t = (b/rho) Lap_xi psi2 = (b/T) e^{tau} Lap_xi psi2
    ok("B7 rho^2 * b Lap u_t = (b/T) e^{tau} Lap_xi psi2: the viscous term grows like e^{tau}")


# ------------------------------------------------------------------ C1, C2
def check_C():
    t, T = sp.symbols('t T', real=True)
    xs = sp.symbols('x1:4', real=True)
    u = sp.Function('u')(t, *xs)
    w = sp.exp(-u)
    lhs = sp.simplify(box(w, t, xs)*sp.exp(u))
    rhs = sp.diff(u, t)**2 - grad2(u, xs) - box(u, t, xs)
    assert sp.simplify(lhs - rhs) == 0
    # null form tilted family
    a = sp.symbols('a1:4', real=True)
    f = T - t + sum(ai*xi for ai, xi in zip(a, xs))
    C = sp.Symbol('C')
    uC = -C*sp.log(f)
    eq = sp.factor(sp.simplify((box(uC, t, xs) - sp.diff(uC, t)**2 + grad2(uC, xs))*f**2))
    sols = sp.solve(eq, C)
    assert 1 in sols, sols
    # non-ODE self-similar blowup
    d = 3
    uQ = -sp.log((T - t)**2 + sum(x**2 for x in xs)/d)
    assert sp.simplify(box(uQ, t, xs) - sp.diff(uQ, t)**2 + grad2(uQ, xs)) == 0
    ok("C1 Nirenberg: Box e^{-u} = e^{-u}(u_t^2-|grad u|^2-Box u); null-form tilt -log(T-t+a.x) "
       "(coefficient 1); -log((T-t)^2+|x|^2/d) exact")
    # model (b) in w = e^{-u}
    lhs = sp.simplify(w*box(w, t, xs) + grad2(w, xs))
    rhs = sp.simplify(sp.exp(-2*u)*(sp.diff(u, t)**2 - box(u, t, xs)))
    assert sp.simplify(lhs - rhs) == 0
    ok("C2 model (b) <=> w Box w + |grad_x w|^2 = 0 for w = e^{-u}")


# ------------------------------------------------------------------ E1
def check_E1():
    X, Y, sp_ = sp.symbols('X Y s', real=True)
    f1 = sp.Function('f1')(X, Y)
    f2 = sp.Function('f2')(X, Y)
    Lam = lambda g: X*sp.diff(g, X) + Y*sp.diff(g, Y)
    Lap = lambda g: sp.diff(g, X, 2) + sp.diff(g, Y, 2)
    gx = lambda g: (sp.diff(g, X), sp.diff(g, Y))
    d = 2
    lhs = (2*sum(a*b for a, b in zip(gx(f1), gx(-Lam(f1) - sp_*f1 + f2)))
           + 2*f2*(Lap(f1) - Lam(f2) - (sp_ + 1)*f2))
    e = sum(a**2 for a in gx(f1)) + f2**2
    flux = (-X*e + 2*f2*gx(f1)[0], -Y*e + 2*f2*gx(f1)[1])
    rhs = (d - 2 - 2*sp_)*e + sp.diff(flux[0], X) + sp.diff(flux[1], Y)
    assert sp.simplify(sp.expand(lhs - rhs)) == 0
    g = sp.Function('g')(X, Y)
    assert sp.simplify(sp.diff(Lam(g), X, 2, Y, 1) - (Lam(sp.diff(g, X, 2, Y, 1)) + 3*sp.diff(g, X, 2, Y, 1))) == 0
    ok("E1 energy identity: dE/dtau = (d-2-2s')E - int_{dB}(|grad f1|^2 + f2^2 - 2 f2 d_r f1) "
       "<= (d-2-2s')E;  D^beta Lam = (Lam+|beta|) D^beta")


# ------------------------------------------------------------------ W1, W2
def check_W():
    t, k, c, b = sp.symbols('t k c b', positive=True)
    xs = sp.symbols('x1:4', real=True)
    p = sp.Function('p')(t, *xs)
    lap = lambda f: sum(sp.diff(f, x, 2) for x in xs)
    west = (1 - 2*k*p)*sp.diff(p, t, 2) - c**2*lap(p) - 2*k*sp.diff(p, t)**2
    # pressure form from (p^2)_tt
    west2 = sp.diff(p, t, 2) - c**2*lap(p) - k*sp.diff(p**2, t, 2)
    assert sp.simplify(west - west2) == 0
    q = 1 - 2*k*p
    qform = (sp.diff(q**2, t, 2) - 2*c**2*lap(q))
    assert sp.simplify(qform + 4*k*west) == 0 or sp.simplify(qform - (-4*k)*west) == 0
    # homogeneous: q^2 = A + B t  => p_t ~ (T-t)^{-1/2}
    A, B, T = sp.symbols('A B T', positive=True)
    qh = sp.sqrt(B*(T - t))
    ph = (1 - qh)/(2*k)
    assert sp.simplify(((1 - 2*k*ph)*sp.diff(ph, t, 2) - 2*k*sp.diff(ph, t)**2)) == 0
    assert sp.simplify(sp.diff(ph, t) - sp.sqrt(B)/(4*k*sp.sqrt(T - t))) == 0
    ok("W1 Westervelt (inviscid) <=> (q^2)_tt = 2c^2 Lap q, q = 1-2kp; homogeneous "
       "q = sqrt(B(T-t)), p -> 1/(2k), p_t = sqrt(B)/(4k sqrt(T-t))")
    # potential forms
    psi = sp.Function('psi')(t, *xs)
    K, beta = sp.symbols('K beta', positive=True)
    grad2p = grad2(psi, xs)
    kuz_rhs = sp.diff(K*sp.diff(psi, t)**2 + grad2p, t)
    west_rhs = sp.diff((K + 1/c**2)*sp.diff(psi, t)**2, t)   # beta/c^2 = (1 + B/2A)/c^2
    diff_ = sp.simplify(kuz_rhs - west_rhs - sp.diff(grad2p - sp.diff(psi, t)**2/c**2, t))
    assert diff_ == 0
    # differentiate Westervelt potential form -> pressure form with p = rho0 psi_t
    rho0 = sp.Symbol('rho0', positive=True)
    Wpot = sp.diff(psi, t, 2) - c**2*lap(psi) - b*lap(sp.diff(psi, t)) - (beta/c**2)*sp.diff(sp.diff(psi, t)**2, t)
    P = rho0*sp.diff(psi, t)
    Wpres = sp.diff(P, t, 2) - c**2*lap(P) - b*lap(sp.diff(P, t)) - (beta/(rho0*c**2))*sp.diff(P**2, t, 2)
    assert sp.simplify(rho0*sp.diff(Wpot, t) - Wpres) == 0
    # Kuznetsov uniform-compression solution
    f = sp.Function('f')(t)
    d = 3
    Bc = sp.Symbol('Bc', real=True)
    psiK = f + Bc*sum(x**2 for x in xs)/(2*d*c**2)
    kuz = (sp.diff(psiK, t, 2) - c**2*lap(psiK) - b*lap(sp.diff(psiK, t))
           - sp.diff(K*sp.diff(psiK, t)**2 + grad2(psiK, xs), t))
    target = sp.diff(sp.diff(f, t) - K*sp.diff(f, t)**2, t) - Bc
    assert sp.simplify(kuz - target) == 0
    # Kuznetsov in r = 1 - 2K psi_t:  (r^2)_tt = 2c^2 Lap r - 8K d_t(grad psi . grad psi_t)
    kuz_inv = (sp.diff(psi, t, 2) - c**2*lap(psi) - sp.diff(K*sp.diff(psi, t)**2 + grad2p, t))
    r_ = 1 - 2*K*sp.diff(psi, t)
    gdot = sum(sp.diff(psi, x)*sp.diff(psi, x, t) for x in xs)
    rform = sp.diff(r_**2, t, 2) - 2*c**2*lap(r_) + 8*K*sp.diff(gdot, t)
    assert sp.simplify(rform + 4*K*sp.diff(kuz_inv, t)) == 0
    ok("W2 Kuznetsov = Westervelt + d_t(|grad psi|^2 - psi_t^2/c^2); d_t(potential Westervelt) "
       "= pressure Westervelt; Kuznetsov uniform compression psi = f(t)+B|x|^2/(2dc^2), "
       "(f'-K f'^2)' = B  (viscous term inactive); (r^2)_tt = 2c^2 Lap r - 8K d_t(grad psi.grad psi_t)")


# ------------------------------------------------------------------ W3, W4
def check_W34():
    # W3 : exact similarity reduction, d = 2, generic Q
    t, x, y, T, c = sp.symbols('t x y T c', positive=True)
    tau, X, Y = sp.symbols('tau X Y', real=True)
    Q = sp.exp(sp.Rational(1, 5)*tau)*(2 + X**2*Y) + sp.cos(X + tau*Y)
    rho = T - t
    phys = {tau: -sp.log(rho), X: x/rho**sp.Rational(3, 4), Y: y/rho**sp.Rational(3, 4)}
    q = sp.sqrt(rho)*Q.subs(phys)
    lhs = sp.diff(q**2, t, 2) - 2*c**2*(sp.diff(q, x, 2) + sp.diff(q, y, 2))
    D = lambda g: X*sp.diff(g, X) + Y*sp.diff(g, Y)
    G = Q**2
    op = lambda g: sp.diff(g, tau) + sp.Rational(3, 4)*D(g)
    red = op(op(G) - G) - 2*c**2*(sp.diff(Q, X, 2) + sp.diff(Q, Y, 2))
    sim = {t: T - sp.exp(-tau), x: X*sp.exp(-sp.Rational(3, 4)*tau), y: Y*sp.exp(-sp.Rational(3, 4)*tau)}
    e = sp.simplify(lhs.subs(sim)*sp.exp(-tau) - red)
    assert e == 0, e
    ok("W3 (q^2)_tt - 2c^2 Lap q = rho^{-1}[(d_tau+3/4 D)(d_tau+3/4 D-1)Q^2 - 2c^2 Lap_eta Q], "
       "q = rho^{1/2} Q(tau, x/rho^{3/4})")
    # W4 : linearisation and eigen-ODE
    lam, r, z = sp.symbols('lambda r z')
    for d in (1, 3):
        for l in (0, 1, 2):
            if d == 1 and l > 1:
                continue
            m1, m2 = sp.Rational(4, 3)*lam, sp.Rational(4, 3)*(lam - 1)
            Hf = sp.Function('H')(r)
            g = r**l*Hf
            Lam = lambda f: r*sp.diff(f, r)
            lapr = lambda f: sp.diff(f, r, 2) + (d - 1)/r*sp.diff(f, r) - l*(l + d - 2)/r**2*f
            ode = sp.simplify(((Lam(Lam(g)) + (m1 + m2)*Lam(g) + m1*m2*g) - lapr(g))/r**l)
            Gz = sp.Function('G')
            e = sp.simplify(ode.subs(Hf, Gz(r**2)).doit().subs(r, sp.sqrt(z)))
            al, be, cc = (l + m1)/2, (l + m2)/2, l + sp.Rational(d, 2)
            G_ = Gz(z)
            hyp = z*(1 - z)*sp.diff(G_, z, 2) + (cc - (al + be + 1)*z)*sp.diff(G_, z) - al*be*G_
            assert sp.simplify(e + 4*hyp) == 0 or sp.simplify(e - 4*hyp) == 0, (d, l, sp.simplify(e/hyp))
    ok("W4a Westervelt modes: (Lam+4lam/3)(Lam+4(lam-1)/3)h = Lap h  <->  2F1 with "
       "alpha=(l+4lam/3)/2, beta=(l+4(lam-1)/3)/2, c=l+d/2")
    # explicit linear modes in physical variables: Z = q0 v solves Z_tt = q0^{-1} Lap Z,
    # q0 = sqrt(A rho); delta(q^2) = 2 Z.
    A = sp.Symbol('A', positive=True)
    rh = sp.Symbol('rho', positive=True)
    xs = sp.symbols('x1:4', real=True)
    d = 3
    q0 = sp.sqrt(A*rh)
    modes = {sp.Integer(1): sp.Integer(1),                 # time translation
             sp.Rational(1, 4): xs[0],                      # tilt of degeneracy set
             sp.Integer(0): rh,                             # amplitude A
             sp.Rational(-1, 2): sum(x**2 for x in xs) + sp.Rational(8, 3)*d/sp.sqrt(A)*rh**sp.Rational(3, 2),
             sp.Rational(-3, 4): rh*xs[0]}
    for lam_, Z in modes.items():
        res = sp.simplify(sp.diff(Z, rh, 2) - sum(sp.diff(Z, x, 2) for x in xs)/q0)
        assert res == 0, (lam_, res)
        # growth relative to background q^2 = A rho in similarity variables:
        # delta(q^2)/rho = rho^{-lam} * function(eta), eta = x/rho^{3/4}
        e_ = sp.Symbol('eta', real=True)
        sub = {xs[0]: e_*rh**sp.Rational(3, 4), xs[1]: e_*rh**sp.Rational(3, 4), xs[2]: 0}
        ratio = sp.simplify((Z/rh).subs(sub)*rh**lam_)
        assert not ratio.has(rh), (lam_, ratio)
    ok("W4b explicit Westervelt modes lam = 1 (time transl.), 1/4 (tilt, l=1), 0 (amplitude), "
       "-1/2 (curvature, l=0), -3/4 (l=1)")


if __name__ == "__main__":
    check_B1()
    check_B2()
    check_B3()
    check_B4()
    check_B5()
    check_B6()
    check_B7()
    check_C()
    check_E1()
    check_W()
    check_W34()
    print("all symbolic checks passed")
