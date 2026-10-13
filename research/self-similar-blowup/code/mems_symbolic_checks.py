#!/usr/bin/env python3
"""Symbolic verification of every formula used in sections/03_mems.tex.

Hyperbolic MEMS (quenching) equation, undamped part:
    w_tt - Delta w = -lam * w^(-q),   q > 1,   sigma = 2/(q+1).
Run:  python3 mems_symbolic_checks.py      (all checks print 'OK' or raise)
"""
import random
import sympy as sp

t, T, lam, q, gam = sp.symbols('t T lambda q gamma', positive=True)
sig = 2 / (q + 1)


def ok(expr, name, samples=6):
    """Check expr == 0: symbolically if possible, otherwise at random admissible
    parameter values (q in (1,6), all other free symbols in (0,1), T = 2) with 30 digits."""
    e = sp.simplify(expr)
    if e == 0:
        print(f"[OK] {name}  (symbolic)")
        return
    free = sorted(e.free_symbols, key=str)
    if e.atoms(sp.Function) - e.atoms(sp.functions.elementary.miscellaneous.Max):
        if any(isinstance(fa, sp.core.function.AppliedUndef) for fa in e.atoms(sp.Function)):
            raise AssertionError(f"{name}: residual {e}")
    rng = random.Random(7)
    for _ in range(samples):
        sub = {}
        for sym in free:
            if sym.name == 'q':
                sub[sym] = sp.Rational(rng.randint(105, 600), 100)
            elif sym.name == 'T':
                sub[sym] = 2
            else:
                sub[sym] = sp.Rational(rng.randint(5, 95), 100)
        v = sp.N(e.subs(sub), 40)
        if abs(v) > 1e-25:
            raise AssertionError(f"{name}: residual {e} = {v} at {sub}")
    print(f"[OK] {name}  (numeric, {samples} random samples)")


# ---------------------------------------------------------------- (a) ODE touchdown
c = (lam / (sig * (1 - sig))) ** (1 / (q + 1))
W = c * (T - t) ** sig
ok(sp.diff(W, t, 2) + lam * W ** (-q), "ODE solution c(T-t)^sigma")
ok(sp.simplify(c ** (q + 1) - lam * (q + 1) ** 2 / (2 * (q - 1))),
   "c^{q+1} = lam (q+1)^2 / (2(q-1))")
# zero energy:  E = W'^2/2 - lam W^{1-q}/(q-1) = 0
E = sp.diff(W, t) ** 2 / 2 - lam * W ** (1 - q) / (q - 1)
ok(sp.simplify(sp.powsimp(sp.expand_power_base(E, force=True), force=True)),
   "ODE solution has zero energy")

# ------------------------------------------------- Lorentz-boosted ODE solution, d=3 general a
x1, x2, x3, a1, a2, a3 = sp.symbols('x1 x2 x3 a1 a2 a3', real=True)
X = [x1, x2, x3]
A = [a1, a2, a3]
adotx = sum(ai * xi for ai, xi in zip(A, X))
a2s = sum(ai ** 2 for ai in A)
Wa = c * (1 - a2s) ** (-sig / 2) * (T - t + adotx) ** sig
res = sp.diff(Wa, t, 2) - sum(sp.diff(Wa, xi, 2) for xi in X) + lam * Wa ** (-q)
# check numerically at random rational points (symbolic simplification is slow)
import random
random.seed(1)
for _ in range(5):
    sub = {q: sp.Rational(random.randint(11, 50), 10), lam: sp.Rational(random.randint(1, 30), 10),
           T: 1, t: sp.Rational(1, 3), x1: sp.Rational(1, 10), x2: sp.Rational(-1, 7), x3: sp.Rational(1, 5),
           a1: sp.Rational(1, 4), a2: sp.Rational(-1, 5), a3: sp.Rational(1, 9)}
    v = sp.N(res.subs(sub), 30)
    assert abs(v) < 1e-20, v
print("[OK] boosted solution c(1-|a|^2)^(-sigma/2) (T-t+a.x)^sigma solves the PDE (d=3, random points)")

# ------------------------------------------- similarity variables and the first-order system
# psi1 = (T-t)^(-sigma) w(t,(T-t) xi),  psi2 = (T-t)^(1-sigma) w_t ; check the static profile eq.
# We verify: if w = (T-t)^sigma psi(x/(T-t)) then w_tt - Delta w + lam w^-q =
#   (T-t)^(sigma-2) [ (delta_ij - xi_i xi_j) d_ij psi - 2(1-sigma) Lambda psi + sigma(1-sigma) psi - lam psi^-q ] * (-1)
xi1, xi2, xi3 = sp.symbols('xi1 xi2 xi3', real=True)
XI = [xi1, xi2, xi3]
s_ = T - t


def profile_op(P):
    """(delta_ij - xi_i xi_j) d_ij P - 2(1-sigma) Lambda P + sigma(1-sigma) P - lam P^-q"""
    LamP = sum(xi * sp.diff(P, xi) for xi in XI)
    return (sum(sp.diff(P, xi, 2) for xi in XI)
            - sum(XI[i] * XI[j] * sp.diff(P, XI[i], XI[j]) for i in range(3) for j in range(3))
            - 2 * (1 - sig) * LamP + sig * (1 - sig) * P - lam * P ** (-q))


# concrete (non-symmetric) test function: the identity is linear in the derivatives of psi,
# so agreement for a generic test function at random points verifies it.
ptest = 2 + sp.sin(xi1) + xi2 ** 2 * xi3 + sp.exp(xi3) / 5 + xi1 * xi2 / 3
wsub = s_ ** sig * ptest.subs({xi1: x1 / s_, xi2: x2 / s_, xi3: x3 / s_}, simultaneous=True)
lhs = sp.diff(wsub, t, 2) - sum(sp.diff(wsub, xi, 2) for xi in X) + lam * wsub ** (-q)
rhs = -s_ ** (sig - 2) * profile_op(ptest).subs({xi1: x1 / s_, xi2: x2 / s_, xi3: x3 / s_}, simultaneous=True)
rng = random.Random(3)
for _ in range(6):
    sub = {q: sp.Rational(rng.randint(105, 600), 100), lam: sp.Rational(rng.randint(1, 30), 10), T: 1,
           t: sp.Rational(rng.randint(1, 50), 100), x1: sp.Rational(rng.randint(-15, 15), 100),
           x2: sp.Rational(rng.randint(-15, 15), 100), x3: sp.Rational(rng.randint(-15, 15), 100)}
    v = sp.N((lhs - rhs).subs(sub), 40)
    assert abs(v) < 1e-25, v
print("[OK] w=(T-t)^sigma psi(x/(T-t)):  w_tt - Delta w + lam w^-q = -(T-t)^(sigma-2) [profile operator](psi)")

# boosted profile solves the profile PDE (symbolic in xi, numeric check)
psia = c * (1 - a2s) ** (-sig / 2) * (1 + sum(ai * xi for ai, xi in zip(A, XI))) ** sig
opa = profile_op(psia)
for _ in range(4):
    sub = {q: sp.Rational(random.randint(11, 50), 10), lam: sp.Rational(random.randint(1, 30), 10),
           xi1: sp.Rational(1, 3), xi2: sp.Rational(-1, 4), xi3: sp.Rational(1, 5),
           a1: sp.Rational(1, 4), a2: sp.Rational(-1, 5), a3: sp.Rational(1, 9)}
    v = sp.N(opa.subs(sub), 30)
    assert abs(v) < 1e-20, v
print("[OK] psi_a(xi) = c(1-|a|^2)^(-sigma/2)(1+a.xi)^sigma solves the profile PDE")

# radial form of the profile ODE
rho, d = sp.symbols('rho d', positive=True)
f = sp.Function('f')(rho)
# Laplacian radial: f'' + (d-1)/rho f'; xi_i xi_j d_ij f = rho^2 f''; Lambda f = rho f'
rad = (1 - rho ** 2) * f.diff(rho, 2) + ((d - 1) / rho - 2 * (1 - sig) * rho) * f.diff(rho) \
    + sig * (1 - sig) * f - lam * f ** (-q)
print("[info] radial profile ODE:", rad)

# ------------------------------------------- first-order system, check against NOTATION.md
# psi1 = (T-t)^s u , psi2 = (T-t)^(s+1) u_t, tau = -log(T-t)+log T, xi = x/(T-t)  (d=1 suffices for the
# structure; the Laplacian part is identical in each variable)
s = sp.symbols('s', real=True)
x = sp.symbols('x', real=True)
u = sp.Function('u')
tau, xi = sp.symbols('tau xi', real=True)
# express t,x in terms of tau, xi (T fixed): T-t = T e^{-tau}, x = xi T e^{-tau}
Psi1 = sp.Function('Psi1')(tau, xi)
Psi2 = sp.Function('Psi2')(tau, xi)
# u(t,x) = (T-t)^{-s} Psi1(tau(t), x/(T-t));  check d_tau Psi1 formula by chain rule
tt = sp.symbols('tt', real=True)
uu = (T - tt) ** (-s) * Psi1.subs({tau: -sp.log(T - tt) + sp.log(T), xi: x / (T - tt)}, simultaneous=True)
# u_t expressed:
ut = sp.diff(uu, tt)
# psi2 := (T-t)^{s+1} u_t  ->  must equal d_tau Psi1 + Lambda Psi1 + s Psi1
expr = sp.simplify((T - tt) ** (s + 1) * ut)
expr = expr.subs({x: xi * (T - tt)}).subs({tt: T - T * sp.exp(-tau)})
expr = sp.simplify(expr)
target = sp.diff(Psi1, tau) + xi * sp.diff(Psi1, xi) + s * Psi1
ok(sp.simplify(expr - target), "psi2 = d_tau psi1 + Lambda psi1 + s psi1  (first equation)")

# --------------------------------------------- linearisation constants
V0 = q * lam * c ** (-q - 1)
ok(V0 - 2 * q * (q - 1) / (q + 1) ** 2, "V0 = q lam c^{-q-1} = 2q(q-1)/(q+1)^2")
aa = sp.symbols('a')
roots = sp.solve(aa ** 2 + aa - V0, aa)
roots = [sp.simplify(r) for r in roots]
print("[info] roots of a^2+a-V0:", roots)
assert any(sp.simplify(r - (1 - sig)) == 0 for r in roots)
assert any(sp.simplify(r - (sig - 2)) == 0 for r in roots)
print("[OK] roots a1 = 1-sigma = (q-1)/(q+1),  a2 = sigma-2 = -2q/(q+1)")

# --------------------------------------------- eigenvalue equation from the first order system
# (mu - L) f = 0 with L f = ( -Lambda f1 - s f1 + f2 , Delta f1 - Lambda f2 - (s+1) f2 + V0 f1 ), s=-sigma
# done in radial-harmonic form: f1 = G(r) Y_l,  d general
r, l, mu = sp.symbols('r l mu')
G = sp.Function('G')(r)
dd = sp.symbols('d')
ss = -sig
Lam_ = lambda F: r * sp.diff(F, r)
Lap_ = lambda F: sp.diff(F, r, 2) + (dd - 1) / r * sp.diff(F, r) - l * (l + dd - 2) / r ** 2 * F
f2 = (mu + ss) * G + Lam_(G)  # from first component
eq2 = mu * f2 - (Lap_(G) - Lam_(f2) - (ss + 1) * f2 + V0 * G)
a_ = mu + ss
radial_eig = (1 - r ** 2) * G.diff(r, 2) + ((dd - 1) / r - 2 * (a_ + 1) * r) * G.diff(r) \
    - l * (l + dd - 2) / r ** 2 * G - (a_ * (a_ + 1) - V0) * G
ok(sp.expand(eq2 + radial_eig), "eigenvalue ODE (1-r^2)G''+((d-1)/r-2(a+1)r)G'-l(l+d-2)/r^2 G-(a(a+1)-V0)G=0, a=mu-sigma")

# hypergeometric reduction G = r^l h(r^2)
z = sp.symbols('z', positive=True)
h = sp.Function('h')
Gs = r ** l * h(r ** 2)
e = radial_eig.subs(G, Gs).doit()
e = sp.simplify(e / r ** l)
e = e.subs(r, sp.sqrt(z))
e = sp.simplify(e)
al = (a_ + l - (1 - sig)) / 2
be = (a_ + l - (sig - 2)) / 2
gm = l + dd / 2
H = sp.Function('H')(z)
hyp = z * (1 - z) * H.diff(z, 2) + (gm - (al + be + 1) * z) * H.diff(z) - al * be * H
# compare e (in terms of h(z)) with 4*hyp
e2 = e.subs(h(z), H).doit()
# replace derivatives
e2 = sp.simplify(e2)
diffexpr = sp.simplify(sp.expand(e2 - 4 * hyp))
ok(diffexpr, "G = r^l h(r^2): h solves 2F1 equation with alpha=(a+l-a1)/2, beta=(a+l-a2)/2, c=l+d/2")
ok(sp.simplify(gm - al - be - ((dd - 1) / 2 + sig - mu)), "c-alpha-beta = (d-1)/2 + sigma - mu")
ok(sp.simplify(al - be - (sig - sp.Rational(3, 2))), "alpha - beta = sigma - 3/2 (never an integer for q>1)")

# eigenvalues
n = sp.symbols('n')
mu1 = sp.solve(sp.Eq(al, -n), mu)[0]
mu2 = sp.solve(sp.Eq(be, -n), mu)[0]
ok(mu1 - (1 - l - 2 * n), "branch alpha=-n: mu = 1 - l - 2n")
ok(sp.simplify(mu2 - (2 * sig - 2 - l - 2 * n)), "branch beta=-n: mu = -(2-2 sigma) - l - 2n")

# -------------------------------------------------- symmetry eigenfunctions
# time translation: d/dT of (T-t)^(-sigma) c (T'-t)^sigma at T'=T  ->  const ; f2=(1-sigma) f1
# Lorentz: d/da_j psi_a at a=0 = c sigma xi_j ; f2 = Lambda f1 - sigma f1 = (1-sigma) xi_j
ok(sp.diff(psia, a1).subs({a1: 0, a2: 0, a3: 0}) - c * sig * xi1, "d/da1 psi_a |_{a=0} = c sigma xi1")

# -------------------------------------------------- Jordan chain checks (inhomogeneous problem)
# (mu-L)g = f  <=>  -[ (1-r^2)G'' + ((d-1)/r-2(a+1)r)G' - l(l+d-2)/r^2 G - (a(a+1)-V0)G ] = F,
# F = f2 + (a+1) f1 + Lambda f1.
F1 = sp.Function('F1')(r)
F2 = sp.Function('F2')(r)
g1 = sp.Function('g1')(r)
g2 = (mu + ss) * g1 + Lam_(g1) - F1
eqi = mu * g2 - (Lap_(g1) - Lam_(g2) - (ss + 1) * g2 + V0 * g1) - F2
rad_g = (1 - r ** 2) * g1.diff(r, 2) + ((dd - 1) / r - 2 * (a_ + 1) * r) * g1.diff(r) \
    - l * (l + dd - 2) / r ** 2 * g1 - (a_ * (a_ + 1) - V0) * g1
ok(sp.expand(eqi + rad_g + F2 + (a_ + 1) * F1 + Lam_(F1)), "inhomogeneous eq: -[radial op] g1 = F2 + (a+1)F1 + Lambda F1")
# mu=1, l=0: f = (1, 1-sigma) -> F = 3-2 sigma
Fval = ((1 - sig) + (1 - sig + 1) * 1).subs({})
ok(sp.simplify(Fval - (3 - 2 * sig)), "mu=1 Jordan source F = 3-2sigma")
# mu=0, l=1: f = (r, (1-sigma) r) -> F = (1-sigma) r + (1-sigma) r + r
ok(sp.simplify(((1 - sig) * r + (-sig + 1) * r + r) - (3 - 2 * sig) * r), "mu=0 Jordan source F = (3-2sigma) r")
# check a(a+1)-V0 at mu=1 vanishes; at mu=0 equals -2(1-sigma)
ok(sp.simplify(((1 - sig) * (2 - sig) - V0)), "a(a+1)-V0 = 0 at mu=1")
ok(sp.simplify((-sig * (1 - sig) - V0) + 2 * (1 - sig)), "a(a+1)-V0 = -2(1-sigma) at mu=0")

# -------------------------------------------------- d=1 conserved energy in rapidity variable
eta = sp.symbols('eta', real=True)
Phi = sp.Function('Phi')(eta)
# w = R^sigma Phi(eta), T-t = R cosh eta, x = R sinh eta.  Check: w_tt - w_xx + lam w^-q = R^{sigma-2}(sigma^2 Phi - Phi'' + lam Phi^-q)
R_ = sp.symbols('R', positive=True)
Phi0 = c * sp.cosh(eta) ** sig
ok(sp.simplify(sp.diff(Phi0, eta, 2) - sig ** 2 * Phi0 - lam * Phi0 ** (-q)).subs({q: 2, lam: 3}),
   "Phi = c cosh(eta)^sigma solves Phi'' = sigma^2 Phi + lam Phi^-q (q=2, lam=3)")
for qq in [sp.Rational(3, 2), 2, 3, 4, sp.Rational(7, 3)]:
    vv = (sp.diff(Phi0, eta, 2) - sig ** 2 * Phi0 - lam * Phi0 ** (-q)).subs({q: qq, lam: sp.Rational(5, 4)})
    for et in [sp.Rational(-3, 2), 0, sp.Rational(7, 10)]:
        assert abs(sp.N(vv.subs(eta, et), 30)) < 1e-25
print("[OK] Phi = c cosh(eta)^sigma solves the rapidity ODE for several q")
# wave operator in Milne coordinates (d=1): verify with explicit chain rule numerically on a test function
Rs, es = sp.symbols('Rs es', positive=True)
tt_, xx_ = sp.symbols('tt_ xx_', real=True)
wtest = (sp.sqrt((T - tt_) ** 2 - xx_ ** 2)) ** sig * sp.exp(sp.Rational(3, 10) * sp.atanh(xx_ / (T - tt_)))
# for Phi = e^{k eta}: Phi'' = k^2 Phi  -> w_tt - w_xx = R^{sigma-2} (sigma^2 - k^2) Phi
lhs1 = sp.diff(wtest, tt_, 2) - sp.diff(wtest, xx_, 2)
Rr = sp.sqrt((T - tt_) ** 2 - xx_ ** 2)
rhs1 = Rr ** (sig - 2) * (sig ** 2 - sp.Rational(9, 100)) * sp.exp(sp.Rational(3, 10) * sp.atanh(xx_ / (T - tt_)))
vv = (lhs1 - rhs1).subs({q: 2, T: 1, tt_: sp.Rational(1, 5), xx_: sp.Rational(1, 3)})
assert abs(sp.N(vv, 30)) < 1e-25
print("[OK] Milne form of the 1D wave operator: w = R^sigma Phi(eta) => w_tt - w_xx = R^{sigma-2}(sigma^2 Phi - Phi'')")

# energy in terms of psi:  2E = (1-xi^2)^{1-sigma} [ -s^2 psi^2 + 2 s xi psi psi' + (1-xi^2) psi'^2 + 2 lam psi^{1-q}/(q-1) ]
X_ = sp.symbols('X_', real=True)
ps = sp.Function('ps')
PhiX = (1 - X_ ** 2) ** (-sig / 2) * ps(X_)
# d/d eta = (1-X^2) d/dX
dPhi = (1 - X_ ** 2) * sp.diff(PhiX, X_)
twoE = dPhi ** 2 - sig ** 2 * PhiX ** 2 + 2 * lam * PhiX ** (1 - q) / (q - 1)
claim = (1 - X_ ** 2) ** (1 - sig) * (-sig ** 2 * ps(X_) ** 2 + 2 * sig * X_ * ps(X_) * ps(X_).diff(X_)
                                      + (1 - X_ ** 2) * ps(X_).diff(X_) ** 2 + 2 * lam * ps(X_) ** (1 - q) / (q - 1))
for qq in [2, 3, sp.Rational(5, 2)]:
    expr_ = (twoE - claim).subs(q, qq)
    expr_ = expr_.subs({ps(X_).diff(X_): sp.Symbol('p1')}).subs({ps(X_): sp.Symbol('p0')})
    for Xv in [sp.Rational(1, 3), sp.Rational(-2, 5)]:
        vv = expr_.subs({X_: Xv, sp.Symbol('p1'): sp.Rational(3, 7), sp.Symbol('p0'): sp.Rational(5, 4), lam: 2})
        assert abs(sp.N(vv, 30)) < 1e-25, vv
print("[OK] energy formula 2E = (1-xi^2)^{1-sigma} Q(xi) in profile variables")

# -------------------------------------------------- free energy identity (top order)
# For the free system with shift s' : dE/dtau = -(s'+1-d/2)(|grad u1|^2+|u2|^2) + boundary <= ... ;
# boundary integrand: -1/2|grad u1|^2 - 1/2 u2^2 + u2 d_r u1 <= -1/2 (d_r u1 - u2)^2
g_, p_, u2_ = sp.symbols('g p u2', real=True)  # g = |tangential grad|^2 >= 0, p = d_r u1
bd = -(p_ ** 2 + g_) / 2 - u2_ ** 2 / 2 + u2_ * p_
ok(sp.expand(bd + (p_ - u2_) ** 2 / 2 + g_ / 2), "boundary term = -(d_r u1 - u2)^2/2 - |grad_T u1|^2/2 <= 0")

# -------------------------------------------------- damped ODE expansion
# W'' + gam W' + lam W^-q = 0,  W = c th^sigma (1 + c1 th + ...), th = T - t
th = sp.symbols('theta', positive=True)
c1 = sp.symbols('c1')
Wd = c * th ** sig * (1 + c1 * th)
# d/dt = -d/dth
expr = sp.diff(Wd, th, 2) - gam * sp.diff(Wd, th) + lam * Wd ** (-q)
expr = sp.series(sp.simplify(expr / (c * th ** (sig - 2))).subs(q, 2), th, 0, 2).removeO()
sol = sp.solve(sp.expand(expr).coeff(th, 1), c1)
print("[info] q=2: c1 =", sol)
assert sp.simplify(sol[0] - gam) == 0
for qq in [sp.Rational(3, 2), 4, 5, sp.Rational(5, 2)]:
    expr = sp.diff(Wd, th, 2) - gam * sp.diff(Wd, th) + lam * Wd ** (-q)
    expr = sp.series(sp.simplify(expr.subs(q, qq) / (c.subs(q, qq) * th ** (sig.subs(q, qq) - 2))), th, 0, 2).removeO()
    sol = sp.solve(sp.expand(expr).coeff(th, 1), c1)
    assert sp.simplify(sol[0] - gam / (3 - qq)) == 0, (qq, sol)
print("[OK] damped ODE: c1 = gamma/(3-q)  (q != 3; q = 3 resonant -> log)")

# general first-order condition (symbolic q): coefficient of th^{sigma-1}
#   c c1 sigma(sigma+1) - gam c sigma - q sigma(1-sigma) c c1 = 0  -> c1 (3-q) = gam
ok(sp.simplify(sig * (sig + 1) - q * sig * (1 - sig) - sig * (3 - q)), "c1*sigma*(3-q) = gamma*sigma identity")

# ODE perturbation exponents: v = th^m, m(m-1) = V0  -> m = 2-sigma or sigma-1
mm = sp.symbols('m')
msol = sp.solve(mm * (mm - 1) - V0, mm)
print("[info] ODE perturbation exponents:", [sp.simplify(x) for x in msol])
assert any(sp.simplify(x - (2 - sig)) == 0 for x in msol) and any(sp.simplify(x - (sig - 1)) == 0 for x in msol)
print("[OK] linearised ODE: (T-t)^{sigma-1} (time translation), (T-t)^{2-sigma} (energy)  -> mu = 1, -(2-2 sigma)")

# Frobenius exponent at rho=1 for the profile ODE: 0 and (d-1)/2 + sigma
yv = sp.symbols('y')
rr = sp.symbols('rr')
ind = 2 * rr * (rr - 1) - ((dd - 1) - 2 * (1 - sig)) * rr
print("[info] indicial roots at rho=1 (profile):", sp.solve(ind, rr))

# damping term in similarity variables: (T-t)^{s+2} * (-gam u_t) = -gam (T-t) psi2
ok(sp.simplify((T - t) ** (-sig + 2) * gam * sp.Symbol('ut') - gam * (T - t) * ((T - t) ** (1 - sig) * sp.Symbol('ut'))),
   "damping contributes -gamma (T-t) psi2 = -gamma T e^{-tau} psi2")

# -------------------------------------------------- singular static solution and its oscillation
Aamp = sp.symbols('A', positive=True)
rr_ = sp.symbols('rr_', positive=True)
psis = Aamp * rr_ ** sig
prof_rad = (1 - rr_ ** 2) * sp.diff(psis, rr_, 2) + ((dd - 1) / rr_ - 2 * (1 - sig) * rr_) * sp.diff(psis, rr_) \
    + sig * (1 - sig) * psis - lam * psis ** (-q)
Asol = (lam / (sig * (sig + dd - 2))) ** (1 / (q + 1))
for qq, dv in [(2, 2), (2, 3), (3, 3), (4, 2), (sp.Rational(5, 2), 3)]:
    vv = prof_rad.subs(Aamp, Asol).subs({q: qq, dd: dv, lam: sp.Rational(7, 5)})
    for rv in [sp.Rational(1, 5), sp.Rational(3, 4)]:
        assert abs(sp.N(vv.subs(rr_, rv), 30)) < 1e-25
print("[OK] psi_* = A rho^sigma, A^{q+1} = lam/(sigma(sigma+d-2)), solves the radial profile ODE (d>=2)")
ok(sp.simplify(q * lam * Asol ** (-q - 1) - q * sig * (sig + dd - 2)), "potential at psi_*: q lam A^{-q-1} = q sigma(sigma+d-2)")
# oscillation (Joseph--Lundgren type) for q=2: q s (s+d-2) > (d-2)^2/4  <=>  2 <= d <= 7
osc = [dv for dv in range(2, 12) if (2 * sp.Rational(2, 3) * (sp.Rational(2, 3) + dv - 2) - sp.Rational((dv - 2) ** 2, 4)) > 0]
print("[info] q=2: oscillatory dimensions", osc)
assert osc == [2, 3, 4, 5, 6, 7]

# -------------------------------------------------- q = 3: Jordan block at mu = -1 (l = 0)
# eigenfunction f_E = (1, sigma-2) = (1,-3/2); chain: radial op g1 = -F, F = f2 + (a+1) f1 + Lambda f1 = -2
q3 = {q: 3}
a_m1 = (-1 - sig).subs(q3)
V03 = (2 * q * (q - 1) / (q + 1) ** 2).subs(q3)
assert sp.simplify(a_m1 * (a_m1 + 1) - V03) == 0
Gc = r ** 2 / dd
radop = (1 - r ** 2) * sp.diff(Gc, r, 2) + ((dd - 1) / r - 2 * (a_m1 + 1) * r) * sp.diff(Gc, r) - (a_m1 * (a_m1 + 1) - V03) * Gc
Fsrc = (sig - 2).subs(q3) + (a_m1 + 1) * 1
ok(sp.simplify(radop + Fsrc), "q=3: g1 = |xi|^2/d solves the Jordan-chain equation at mu=-1 (eigenvalue -1 = -omega_0 not semisimple)")

# -------------------------------------------------- Jordan chain at mu=1 (l=0) and mu=0 (l=1): first-order reduction
# P solves (1-r^2)P' + ((D-1)/r - 2(2-sigma) r) P = -(3-2 sigma); D = d (mu=1, P = G'), D = d+2 (mu=0, G = rK, P = K')
Dd = sp.symbols('D', positive=True)
Pf = sp.Function('P')(r)
eqP = (1 - r ** 2) * Pf.diff(r) + ((Dd - 1) / r - 2 * (2 - sig) * r) * Pf
Ifac = r ** (Dd - 1) * (1 - r ** 2) ** ((2 - sig) - (Dd - 1) / 2)
# d/dr (I P) = I/(1-r^2) * [ (1-r^2) P' + (...) P ]
lhsI = sp.diff(Ifac * Pf, r)
ok(sp.simplify(lhsI - Ifac / (1 - r ** 2) * eqP), "integrating factor r^{D-1}(1-r^2)^{(2-sigma)-(D-1)/2} for the Jordan-chain ODE")
G1 = sp.Function('G1')(r)
mu1op = (1 - r ** 2) * G1.diff(r, 2) + ((dd - 1) / r - 2 * (2 - sig) * r) * G1.diff(r)
ok(sp.simplify(mu1op - eqP.subs(Pf, G1.diff(r)).subs(Dd, dd).doit()), "mu=1, l=0: radial op = first-order op on G' with D=d")
print("all additional checks passed")
print("all checks passed")
