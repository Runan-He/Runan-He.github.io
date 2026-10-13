#!/usr/bin/env python3
"""therm_symbolic.py -- symbolic checks for Section 4 (hyperbolic thermal explosion).

Every identity used in sections/04_thermal.tex that is purely algebraic is
verified here with sympy.  Run:  python3 code/therm_symbolic.py
Each check prints OK or raises AssertionError.

Contents
  1. Cattaneo derivation: (1+tau d_t)(energy balance) identity, FK variable.
  2. Reduced model u_tt + g u_t - Lap u = delta e^u: similarity system (d=2,
     generic u), ODE blowup, Lorentz-boosted family, scaling.
  3. Hypergeometric reduction of
        (delta_ij - xi_i xi_j) d_ij g - B Lambda g - C g = 0,  g = r^l Y_l h(r^2)
     for general d, l (radial ODE form), and the alpha, beta of both models.
  4. Jordan-chain equations (reduced, lambda=1,0; full, lambda=-1) and the
     explicit generalised eigenfunction of the full model.
  5. 1D Liouville equation: general solution, Schwarzian invariant, local
     expansion at a blowup point -> boosted profile + O(T-t) (lambda=-1 mode).
  6. Full Cattaneo model: ODE reduction (1+tau d_t)(theta' - delta e^theta)=0,
     normal form (d_t+g)(u_t - e^u) = Lap u, similarity system, tilted family,
     Fuchsian log coefficient.
  7. Damped reduced ODE: resonant log term in the blowup expansion.
"""
import sympy as sp

ok = lambda msg: print(f"[OK] {msg}")


def is_zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


# ---------------------------------------------------------------------------
# 1. Cattaneo derivation
# ---------------------------------------------------------------------------
t, x, y, z = sp.symbols('t x y z', real=True)
tau_r = sp.symbols('tau_r', positive=True)
theta = sp.Function('theta')(t, x)
q = sp.Function('q')(t, x)
f = sp.Function('f')
# 1D: energy balance  theta_t + q_x = f(theta);  Cattaneo  q + tau q_t = -theta_x
E = sp.diff(theta, t) + sp.diff(q, x) - f(theta)
C = q + tau_r * sp.diff(q, t) + sp.diff(theta, x)
# (1+tau d_t)E - d_x C  should eliminate q entirely:
expr = E + tau_r * sp.diff(E, t) - sp.diff(C, x)
target = (tau_r * sp.diff(theta, t, 2) + sp.diff(theta, t) - sp.diff(theta, x, 2)
          - f(theta) - tau_r * sp.diff(f(theta), t))
assert is_zero(expr - target)
ok("Cattaneo: (1+tau d_t)(balance) - div(flux law) = full hyperbolic model")

# Frank-Kamenetskii variable: exp(-E/RT) with T = T0(1+eps*theta), eps = R T0/E
eps, th = sp.symbols('epsilon theta', positive=True)
Ea_over_RT0 = 1 / eps
lhs = sp.exp(-Ea_over_RT0 / (1 + eps * th))
rhs = sp.exp(-Ea_over_RT0) * sp.exp(th / (1 + eps * th))
assert is_zero(sp.log(lhs) - sp.log(rhs))
ok("FK variable: exp(-E/RT) = exp(-E/RT0) exp(theta/(1+eps theta)) exactly")

# ---------------------------------------------------------------------------
# 2. Reduced model: similarity system in d = 2 for a generic function u
# ---------------------------------------------------------------------------
T, delta, gam = sp.symbols('T delta gamma', positive=True)
u = sp.Function('u')(t, x, y)
s_ = T - t
psi1 = u + 2 * sp.log(s_)
psi2 = s_ * sp.diff(u, t)
Dtau = lambda F: s_ * sp.diff(F, t) - x * sp.diff(F, x) - y * sp.diff(F, y)
Lam = lambda F: x * sp.diff(F, x) + y * sp.diff(F, y)
Lapxi = lambda F: s_**2 * (sp.diff(F, x, 2) + sp.diff(F, y, 2))
r1 = Dtau(psi1) - (-Lam(psi1) + psi2 - 2)
assert is_zero(r1)
pde = sp.diff(u, t, 2) + gam * sp.diff(u, t) - sp.diff(u, x, 2) - sp.diff(u, y, 2) - delta * sp.exp(u)
r2 = Dtau(psi2) - (Lapxi(psi1) - Lam(psi2) - psi2 + delta * sp.exp(psi1) - gam * s_ * psi2)
assert is_zero(r2 - s_**2 * pde)
ok("reduced model: similarity system (incl. damping -gamma (T-t) psi2) verified, d=2")

# ODE blowup and boosted family in d = 3
a1, a2, a3 = sp.symbols('a1 a2 a3', real=True)
A0 = sp.sqrt(1 + a1**2 + a2**2 + a3**2)
cc = sp.log(2 / delta)
ub = cc - 2 * sp.log(A0 * (T - t) + a1 * x + a2 * y + a3 * z)
res = (sp.diff(ub, t, 2) - sp.diff(ub, x, 2) - sp.diff(ub, y, 2) - sp.diff(ub, z, 2)
       - delta * sp.exp(ub))
assert is_zero(res)
ok("boosted ODE blowup u_{T,a} = log(2/delta) - 2 log(A0(T-t) + a.x), A0^2 = 1+|a|^2, solves u_tt - Lap u = delta e^u (d=3)")
# its similarity representation is static
xi1, xi2, xi3, tau = sp.symbols('xi1 xi2 xi3 tau', real=True)
sub = {x: (T - t) * xi1, y: (T - t) * xi2, z: (T - t) * xi3}
P1 = sp.simplify((ub + 2 * sp.log(T - t)).subs(sub))
P2 = sp.simplify(((T - t) * sp.diff(ub, t)).subs(sub))
print("   Psi_a,1 =", sp.simplify(sp.expand_log(P1, force=True)), "  Psi_a,2 =", P2)
assert sp.diff(sp.expand_log(P1, force=True), t) == 0 and sp.simplify(sp.diff(P2, t)) == 0
ok("boosted profiles are tau-independent in similarity variables")
# scaling symmetry
mu = sp.symbols('mu', positive=True)
U = sp.Function('U')
v = U(mu * t, mu * x) + 2 * sp.log(mu)
lhs_s = sp.diff(v, t, 2) - sp.diff(v, x, 2) - delta * sp.exp(v)
X1, X2 = sp.symbols('X1 X2')
Ufun = sp.Function('U')(X1, X2)
lhs_U = (sp.diff(Ufun, X1, 2) - sp.diff(Ufun, X2, 2) - delta * sp.exp(Ufun)).subs({X1: mu * t, X2: mu * x})
assert is_zero(lhs_s.doit() - mu**2 * lhs_U.doit())
ok("scaling u -> u(mu t, mu x) + 2 log mu is a symmetry (and gamma -> mu gamma with damping)")

# ---------------------------------------------------------------------------
# 3. Hypergeometric reduction
# ---------------------------------------------------------------------------
r, rho = sp.symbols('r rho', positive=True)
d, l, B, Cc = sp.symbols('d ell B C')
h = sp.Function('h')
fr = r**l * h(r**2)
# operator on f(r) Y_l:  (1-r^2) f'' + (d-1)/r f' - l(l+d-2)/r^2 f - B r f' - C f
op = ((1 - r**2) * sp.diff(fr, r, 2) + (d - 1) / r * sp.diff(fr, r)
      - l * (l + d - 2) / r**2 * fr - B * r * sp.diff(fr, r) - Cc * fr)
cpar = l + d / 2
sab = l + (B - 1) / 2          # alpha + beta
pab = (l * (l - 1) + B * l + Cc) / 4  # alpha * beta
H = sp.Function('H')
hyp = (rho * (1 - rho) * sp.diff(H(rho), rho, 2) + (cpar - (sab + 1) * rho) * sp.diff(H(rho), rho)
       - pab * H(rho))
hyp_r = 4 * r**l * hyp.subs(rho, r**2)
# compare after replacing derivatives of h(r^2)
w = sp.symbols('w', positive=True)
op2 = op.doit()
diffexpr = sp.simplify((op2 - hyp_r.doit().subs(H(r**2), h(r**2))).rewrite(sp.Derivative))
# robust check: substitute a concrete test function h = exp(k rho) etc.
for test in [sp.exp(sp.Rational(1, 3) * w), w**3 + 2 * w + 1, sp.sin(w)]:
    ftest = r**l * test.subs(w, r**2)
    lhs = ((1 - r**2) * sp.diff(ftest, r, 2) + (d - 1) / r * sp.diff(ftest, r)
           - l * (l + d - 2) / r**2 * ftest - B * r * sp.diff(ftest, r) - Cc * ftest)
    Ht = test.subs(w, rho)
    rhs = 4 * r**l * (rho * (1 - rho) * sp.diff(Ht, rho, 2) + (cpar - (sab + 1) * rho) * sp.diff(Ht, rho)
                      - pab * Ht).subs(rho, r**2)
    assert is_zero(lhs - rhs), test
ok("ell-sector reduces to 4 r^l [rho(1-rho)h'' + (c-(a+b+1)rho)h' - ab h], c=l+d/2, a+b=l+(B-1)/2, ab=(l(l-1)+Bl+C)/4")

lam = sp.symbols('lambda')
# reduced model: B = 2(lam+1), C = lam(lam+1) - 2
al_R, be_R = (lam + l - 1) / 2, (lam + l + 2) / 2
assert is_zero(al_R + be_R - sab.subs(B, 2 * (lam + 1)))
assert is_zero(al_R * be_R - pab.subs({B: 2 * (lam + 1), Cc: lam * (lam + 1) - 2}))
ok("reduced: alpha=(lam+l-1)/2, beta=(lam+l+2)/2  =>  lam = 1-l-2n  or  -2-l-2n")
# full model: (Lambda+lam)^2 g - Lap g - g = 0  <=>  B = 2 lam + 1, C = lam^2 - 1
al_F, be_F = (lam + l - 1) / 2, (lam + l + 1) / 2
assert is_zero(al_F + be_F - sab.subs(B, 2 * lam + 1))
assert is_zero(al_F * be_F - pab.subs({B: 2 * lam + 1, Cc: lam**2 - 1}))
ok("full:    alpha=(lam+l-1)/2, beta=(lam+l+1)/2  =>  lam = 1-l-2n  or  -1-l-2n")
# c - alpha - beta
print("   reduced c-a-b =", sp.simplify(cpar - al_R - be_R), "; full c-a-b =", sp.simplify(cpar - al_F - be_F))

# pencils: check (Lambda+lam)(Lambda+lam+1)g - Lap g - 2g corresponds to B,C above (d=2 generic g)
g = sp.Function('g')(x, y)
LamG = lambda F: x * sp.diff(F, x) + y * sp.diff(F, y)
Lap = lambda F: sp.diff(F, x, 2) + sp.diff(F, y, 2)
pencilR = LamG(LamG(g) + (lam + 1) * g) + lam * (LamG(g) + (lam + 1) * g) - Lap(g) - 2 * g
formR = ((1 - x**2) * sp.diff(g, x, 2) - 2 * x * y * sp.diff(g, x, y) + (1 - y**2) * sp.diff(g, y, 2)
         - 2 * (lam + 1) * LamG(g) - (lam * (lam + 1) - 2) * g)
assert is_zero(pencilR + formR)
pencilF = LamG(LamG(g) + lam * g) + lam * (LamG(g) + lam * g) - Lap(g) - g
formF = ((1 - x**2) * sp.diff(g, x, 2) - 2 * x * y * sp.diff(g, x, y) + (1 - y**2) * sp.diff(g, y, 2)
         - (2 * lam + 1) * LamG(g) - (lam**2 - 1) * g)
assert is_zero(pencilF + formF)
ok("quadratic pencils of both linearised first-order systems match the (B,C) forms")

# ---------------------------------------------------------------------------
# 4. Jordan chains
# ---------------------------------------------------------------------------
# integrating factor identity used in the semisimplicity lemma (alpha = 0):
c_, b_ = sp.symbols('c beta', positive=True)
Hh = sp.Function('Hh')(rho)
lhs = sp.diff(rho**c_ * (1 - rho)**(b_ + 1 - c_) * sp.diff(Hh, rho), rho)
bracket = (rho * (1 - rho) * sp.diff(Hh, rho, 2) + (c_ - (b_ + 1) * rho) * sp.diff(Hh, rho))
quot = sp.expand(sp.powsimp(sp.expand(lhs / (rho**(c_ - 1) * (1 - rho)**(b_ - c_))), force=True))
assert is_zero(quot - sp.expand(bracket))
ok("d/drho[rho^c (1-rho)^(beta+1-c) h'] = rho^(c-1)(1-rho)^(beta-c) * (hypergeometric op with alpha=0)")

# Full model: explicit generalised eigenfunction at lambda=-1:
#   v = -(T-t) log(T-t) - |x|^2/(d (T-t)) solves v_tt - Lap v = d_t( v/(T-t) )
s = sp.symbols('s', positive=True)
for dim in [1, 2, 3, 5]:
    X = sp.symbols('x1:%d' % (dim + 1), real=True)
    r2 = sum(xx**2 for xx in X)
    vv = -(T - t) * sp.log(T - t) - r2 / (dim * (T - t))
    res = sp.diff(vv, t, 2) - sum(sp.diff(vv, xx, 2) for xx in X) - sp.diff(vv / (T - t), t)
    assert is_zero(res)
ok("full model: Jordan chain at lambda=-1: v = -(T-t)log(T-t) - |x|^2/(d(T-t)) solves the linearised eq (d=1,2,3,5)")
# and v = T - t is the eigenfunction (lambda = -1, g = 1)
vv = T - t
assert is_zero(sp.diff(vv, t, 2) - sp.diff(vv / (T - t), t))
ok("full model: eigenfunction lambda=-1, g=1 (v = T-t)")

# Reduced model: lambda=-1 eigenfunction in l=0: h = 2F1(-1, 1/2 ; d/2 ; rho) = 1 - rho/(d/2)*(1/2)
for dim in [1, 2, 3, 4]:
    X = sp.symbols('x1:%d' % (dim + 1), real=True)
    r2 = sum(xx**2 for xx in X)
    gfun = 1 - r2 / dim     # 2F1(-1,1/2;d/2;rho) = 1 - rho/d
    vv = (T - t) * gfun.subs({xx: xx / (T - t) for xx in X})
    res = sp.diff(vv, t, 2) - sum(sp.diff(vv, xx, 2) for xx in X) - 2 * vv / (T - t)**2
    assert is_zero(res)
ok("reduced model: lambda=-1, l=0 eigenfunction g = 1 - |xi|^2/d (v=(T-t)g) verified, d=1..4")

# ---------------------------------------------------------------------------
# 5. Liouville equation in 1D
# ---------------------------------------------------------------------------
Fq, Gq = sp.Function('F'), sp.Function('G')
sig, eta = t + x, t - x
uL = sp.log(8 / delta * sp.diff(Fq(sig), x) * (-sp.diff(Gq(eta), x)) / (Fq(sig) + Gq(eta))**2)
# note: d/dx F(t+x) = F'(sig), -d/dx G(t-x) = G'(eta)
res = sp.diff(uL, t, 2) - sp.diff(uL, x, 2) - delta * sp.exp(uL)
assert is_zero(sp.simplify(res))
ok("Liouville: u = log[ 8 F'(t+x) G'(t-x) / (delta (F(t+x)+G(t-x))^2) ] solves u_tt - u_xx = delta e^u (F'G'>0)")
# Schwarzian invariant: u_ss - u_s^2/2 = S(F) with s = sigma
S1, E1 = sp.symbols('S1 E1')
Fs, Ge = sp.Function('Fs')(S1), sp.Function('Ge')(E1)
uSE = sp.log(8 / delta * sp.diff(Fs, S1) * sp.diff(Ge, E1) / (Fs + Ge)**2)
inv = sp.diff(uSE, S1, 2) - sp.diff(uSE, S1)**2 / 2
schw = sp.diff(Fs, S1, 3) / sp.diff(Fs, S1) - sp.Rational(3, 2) * (sp.diff(Fs, S1, 2) / sp.diff(Fs, S1))**2
assert is_zero(sp.simplify(inv - schw))
assert is_zero(sp.simplify(4 * sp.diff(uSE, S1, E1) - delta * sp.exp(uSE)))
ok("Liouville: u_{sigma sigma} - u_sigma^2/2 = Schwarzian S(F)(sigma)  (null coords, 4u_{sigma eta} = delta e^u)")

# local expansion at a blowup point (T,0): F(T)+G(T)=0, F'(T)=p, G'(T)=q, pq>0
p, qq, F2, G2, F3, G3, F0 = sp.symbols('p q F2 G2 F3 G3 F0', real=True)
xi = sp.symbols('xi', real=True)
sig_ = s * (xi - 1)       # sigma - T
eta_ = -s * (1 + xi)      # eta - T
Fser = F0 + p * sig_ + F2 * sig_**2 / 2 + F3 * sig_**3 / 6
Gser = -F0 + qq * eta_ + G2 * eta_**2 / 2 + G3 * eta_**3 / 6
dF = p + F2 * sig_ + F3 * sig_**2 / 2
dG = qq + G2 * eta_ + G3 * eta_**2 / 2
Phi_over_s = sp.expand((Fser + Gser) / s)
psi1_L = sp.log(8 / delta) + sp.log(dF * dG) - 2 * sp.log(Phi_over_s)
ser = sp.series(psi1_L, s, 0, 2).removeO()
P0 = sp.simplify(ser.subs(s, 0))
vboost = -(p - qq) / (p + qq)     # a/A0
A0v = 1 / sp.sqrt(1 - vboost**2)
Psi_a1 = sp.log(2 / delta) - 2 * sp.log(A0v * (1 + vboost * xi))
chk = sp.simplify(sp.exp(P0 - Psi_a1))
assert sp.simplify(chk - 1) == 0, chk
ok("Liouville: leading term at a blowup point = boosted ODE profile with a/A0 = -(F'-G')/(F'+G') = slope of blowup curve")
g1 = sp.simplify(sp.diff(ser, s))
print("   O(s) correction g1(xi) =", sp.factor(g1))
# symmetric case p = q: g1 should be proportional to the reduced lambda=-1 eigenfunction 1 - xi^2
g1sym = sp.factor(sp.simplify(g1.subs(qq, p)))
print("   p=q:  g1 =", g1sym)
# check (L_a + 1) annihilates the O(s) term of (psi1, psi2): construct psi2 expansion too
uLs = psi1_L - 2 * sp.log(s)                       # u in terms of (s, xi), with x = s*xi
# psi2 = (T-t) u_t ; in (s,xi): u_t = -u_s|_x, and at fixed x: xi = x/s
Xs = sp.symbols('Xs', real=True)
u_sx = uLs.subs(xi, Xs / s)
psi2_L = sp.simplify((s * (-sp.diff(u_sx, s))).subs(Xs, s * xi))
ser2 = sp.series(psi2_L, s, 0, 2).removeO()
g2 = sp.simplify(sp.diff(ser2, s))
Q0 = sp.simplify(ser2.subs(s, 0))
# linearised operator around the boosted profile, eigenvalue -1:
#   -g1 = -xi g1' + g2 ;   -g2 = g1'' - xi g2' - g2 + delta e^{Psi_a1} g1
e1 = sp.simplify(-g1 - (-xi * sp.diff(g1, xi) + g2))
e2 = sp.simplify(-g2 - (sp.diff(g1, xi, 2) - xi * sp.diff(g2, xi) - g2 + delta * sp.exp(P0) * g1))
assert e1 == 0 and e2 == 0, (e1, e2)
ok("Liouville: the O(T-t) correction is an eigenfunction of L_a with eigenvalue -1 (exact convergence rate e^{-tau})")
# identify g1 with the curvature of the blowup curve F(T(x)+x) + G(T(x)-x) = 0 at x = 0
Tp = (qq - p) / (p + qq)                                   # T'(0)
Tpp = -(F2 * (1 + Tp)**2 + G2 * (1 - Tp)**2) / (p + qq)   # T''(0) (implicit differentiation)
A0sq = 1 / (1 - Tp**2)
assert sp.simplify(Tp - vboost) == 0
assert sp.simplify(g1 - Tpp * A0sq * (1 - xi**2) / (1 + Tp * xi)) == 0
ok("Liouville: g1 = T''(x0) A0^2 (1-xi^2)/(1+T'(x0) xi): the e^{-tau} correction is driven by the curvature of the blowup curve")

# ---------------------------------------------------------------------------
# 6. Full Cattaneo model
# ---------------------------------------------------------------------------
tht = sp.Function('theta')(t)
expr = tau_r * sp.diff(tht, t, 2) + sp.diff(tht, t) - delta * sp.exp(tht) - tau_r * delta * sp.exp(tht) * sp.diff(tht, t)
expr2 = (sp.diff(tht, t) - delta * sp.exp(tht)) + tau_r * sp.diff(sp.diff(tht, t) - delta * sp.exp(tht), t)
assert is_zero(expr - expr2)
thsol = -sp.log(delta * (T - t))
assert is_zero(expr.subs(tht, thsol).doit())
ok("full ODE: tau th'' + th' - delta e^th - tau delta e^th th' = (1+tau d_t)(th' - delta e^th); th = -log(delta(T-t)) exact")
# -2 log ansatz is inconsistent: leading power mismatch
cst = sp.symbols('c0', real=True)
th2 = cst - 2 * sp.log(T - t)
lead = sp.simplify(expr.subs(tht, th2).doit() * (T - t)**3)
print("   (T-t)^3 * residual of c-2log(T-t) ansatz at t->T:", sp.simplify(sp.limit(lead, t, T)))
# normal form: t = sqrt(tau) t', theta = u - log(sqrt(tau) delta)
gam_ = 1 / sp.sqrt(tau_r)
uN = sp.Function('uN')(t, x)
lhsN = (sp.diff(sp.diff(uN, t) - sp.exp(uN), t) + gam_ * (sp.diff(uN, t) - sp.exp(uN)) - sp.diff(uN, x, 2))
expanded = sp.diff(uN, t, 2) + gam_ * sp.diff(uN, t) - sp.diff(uN, x, 2) - sp.exp(uN) * sp.diff(uN, t) - gam_ * sp.exp(uN)
assert is_zero(lhsN - expanded)
# check that theta(t,x) = uN(t/sqrt(tau), x) - log(sqrt(tau) delta) maps the full model to the normal form
Uf = sp.Function('Uf')
thN = Uf(t / sp.sqrt(tau_r), x) - sp.log(sp.sqrt(tau_r) * delta)
fullpde = (tau_r * sp.diff(thN, t, 2) + sp.diff(thN, t) - sp.diff(thN, x, 2) - delta * sp.exp(thN)
           - tau_r * delta * sp.exp(thN) * sp.diff(thN, t))
Tt = sp.symbols('Tt')
Uv = sp.Function('Uf')(Tt, x)
normal = (sp.diff(Uv, Tt, 2) + gam_ * sp.diff(Uv, Tt) - sp.diff(Uv, x, 2) - sp.exp(Uv) * sp.diff(Uv, Tt)
          - gam_ * sp.exp(Uv)).subs(Tt, t / sp.sqrt(tau_r))
assert is_zero(sp.simplify(fullpde.doit() - normal.doit()))
ok("full model in normal form: (d_t + gamma)(u_t - e^u) = Lap u, gamma = tau_r^{-1/2} (after t -> sqrt(tau_r) t, shift)")
# similarity system (d = 2, generic u):  psi1 = u + log(T-t), psi2 = (T-t) u_t
psi1F = u + sp.log(s_)
psi2F = s_ * sp.diff(u, t)
r1 = Dtau(psi1F) - (-Lam(psi1F) + psi2F - 1)
assert is_zero(r1)
pdeF = (sp.diff(u, t, 2) - sp.diff(u, x, 2) - sp.diff(u, y, 2) - sp.exp(u) * sp.diff(u, t)
        + gam * (sp.diff(u, t) - sp.exp(u)))
r2 = Dtau(psi2F) - (Lapxi(psi1F) - Lam(psi2F) - psi2F + sp.exp(psi1F) * psi2F
                    - gam * s_ * (psi2F - sp.exp(psi1F)))
assert is_zero(r2 - s_**2 * pdeF)
ok("full model: similarity system psi1 = u+log(T-t), psi2=(T-t)u_t verified (with gamma term)")
# tilted family (gamma = 0)
v1, v2, v3 = sp.symbols('v1 v2 v3', real=True)
ut = sp.log(1 - v1**2 - v2**2 - v3**2) - sp.log(T - t + v1 * x + v2 * y + v3 * z)
res = (sp.diff(ut, t, 2) - sp.diff(ut, x, 2) - sp.diff(ut, y, 2) - sp.diff(ut, z, 2)
       - sp.exp(ut) * sp.diff(ut, t))
assert is_zero(res)
ok("full model (gamma=0): tilted family u_v = log(1-|v|^2) - log(T-t+v.x) is exact (d=3)")
resg = sp.simplify(sp.diff(ut, t) - sp.exp(ut))
print("   gamma-defect u_t - e^u on tilted family:", sp.factor(resg))
# Fuchsian log coefficient: u = -log(phi) + w0 + a phi log phi, phi = T(x) - t (1D)
Tx = sp.Function('Tx')(x)
phi = Tx - t
w0 = sp.log(1 - sp.diff(Tx, x)**2)
acoef = sp.symbols('acoef')
uF = -sp.log(phi) + w0 + acoef * phi * sp.log(phi)
resF = sp.diff(uF, t, 2) - sp.diff(uF, x, 2) - sp.exp(uF) * sp.diff(uF, t)
# leading orders in phi: substitute t = Tx - phi and expand in phi
P = sp.symbols('P', positive=True)
resP = sp.expand(sp.simplify(resF.subs(t, Tx - P)))
lead1 = sp.limit(sp.simplify(resP * P), P, 0)
sol = sp.solve(sp.Eq(sp.simplify(lead1), 0), acoef)
print("   Fuchsian log coefficient a(x) =", sp.simplify(sol[0]))
assert is_zero(sol[0] + sp.diff(Tx, x, 2) / (2 * (1 - sp.diff(Tx, x)**2)))
ok("full model: curved blowup surface forces (T(x)-t)log(T(x)-t) terms, a = -T''/(2(1-T'^2))  (Jordan block at -1)")

# ---------------------------------------------------------------------------
# 7. Damped reduced ODE: u_tt + gamma u_t = delta e^u, s = T - t
# ---------------------------------------------------------------------------
b1, b2, b3 = sp.symbols('b1 b2 b3')
wS = b1 * s + s**2 * (b2 + b3 * sp.log(s))
uS = sp.log(2 / delta) - 2 * sp.log(s) + wS
# d/dt = -d/ds
resS = sp.diff(uS, s, 2) - gam * sp.diff(uS, s) - delta * sp.exp(uS)
resS = sp.expand(sp.series(sp.expand(resS * s**2), s, 0, 3).removeO())
c1 = resS.coeff(s, 1)
sol1 = sp.solve(c1, b1)
assert sol1 == [gam]
c2 = sp.expand(resS.coeff(s, 2).subs(b1, gam))
c2 = sp.simplify(c2)
sol3 = sp.solve(sp.Eq(c2.coeff(sp.log(s), 0), 0), b3)
print("   damped ODE: b1 =", sol1[0], ", b3 =", sol3)
assert is_zero(sol3[0] - 2 * gam**2 / 3)
ok("damped reduced ODE: u = log(2/delta) - 2log s + gamma s + (2 gamma^2/3) s^2 log s + b2 s^2 + ...")

print("\nAll symbolic checks passed.")
