# The single layer of LM26 = Luna-Velasco & Mengual, arXiv:2610.34504, in the variables of the manuscript.
#
# Notation of the manuscript: rho_m intermediate density, mu_m = 1 + A rho_m, c = (1-A^2) kappa (opening of LM26),
# sigma = sigma_0 and gamma = gamma_0 per unit arclength, frame (tau_0, n_0) with n_0 = tau_0^perp,
# perp(x1,x2) = (-x2,x1), nabla^perp = (-d_2, d_1).  Quantities of a state (rho,u,m), eq. (3.5) of the manuscript:
#     mu = 1 + A rho,  v = u + A m,  d = 2[(1-A^2) m - (rho-A) v]/(1-rho^2),  D = d + 2 e2,  a = e2 + (A/mu) v.
# LM26 references: (4.6) momentum ansatz, (4.10) Darcy law in the mixing zone, (4.14) exact hull condition,
# (7.9)-(7.11) leading-order state, (7.12) and the proof of Prop. 7.4 (normal part), (3.7b) left relation (tangential
# part). Equation numbers of the manuscript refer to its revision of 9 Nov 2026.
#
# (a) exact identities in LM26's mixing zone (theta = A);
# (b) the leading-order state of LM26's single layer at t = 0: its tilt delta_LM, its exact hull range, and the
#     identity D.(D-2a) = (2(1-A^2)/mu_m) Phi(kappa) + (delta - delta_*)^2 of Lemma 5.2;
# (c) numerics on the overturned example z_0 = (alpha - 1.5 sin alpha, 0.8 sin alpha);
# (d) the relation between the flux mismatch M of the manuscript and LM26's d_alpha G_+ - d_alpha G_-.
import sympy as sp
import numpy as np

def perp(w): return sp.Matrix([-w[1], w[0]])
def dot(x, y): return (x.T*y)[0]
def zero(expr): return sp.simplify(sp.expand(expr)) == 0

# ---------------------------------------------------------------- (a) exact identities in LM26's mixing zone
A = sp.symbols('A', real=True)
u1, u2, g1, g2, p1, p2 = sp.symbols('u1 u2 g1 g2 p1 p2', real=True)
e2 = sp.Matrix([0, 1])
u = sp.Matrix([u1, u2]); gradg = sp.Matrix([g1, g2]); gperp = perp(gradg)
w = A*u + e2
m = A*u - gperp - w/2                                  # LM26 (4.6) in the mixing zone, theta = A
gradp = -(1 + A**2/2)*u + A*gperp - (A/2)*e2           # LM26 (4.10), theta = A, C(A) = A/2, solved for grad p
v = u + A*m
assert all(zero(x) for x in v - (-gradp - A*e2))       # Darcy field v = u + A m = -grad p - A e2
rho = A; mu = 1 + A*rho
D = 2*((1 - A**2)*m - (rho - A)*v)/(1 - rho**2) + 2*e2
a = e2 + (A/mu)*v
assert all(zero(x) for x in D - (w - 2*gperp))                     # D = (Au+e2) - 2 grad^perp g
assert all(zero(x) for x in a - (e2 - A*gradp)/(1 + A**2))         # a = (e2 - A grad p)/(1+A^2)
assert zero(dot(D, D - 2*a) - (4*dot(gradg, gradg) - dot(w, w))/(1 + A**2))
# LM26 Sec. 4.7: 2(1-A theta)(m - theta u) + (1-theta^2)(Au+e2) = -2(1-A^2) grad^perp g, so (O) <=> (4.14)
assert all(zero(x) for x in 2*(1 - A*A)*(m - A*u) + (1 - A**2)*w + 2*(1 - A**2)*gperp)
print("(a) in LM26's mixing zone (theta = A):")
print("    v = u + A m = -grad p - A e2;  D = (Au+e2) - 2 grad^perp g;  a = (e2 - A grad p)/(1+A^2);")
print("    D.(D-2a) = (4|grad g|^2 - |Au+e2|^2)/(1+A^2): LM26's exact hull 2|grad g| < |Au+e2| (4.14) is |D-a| < |a|")

# ---------------------------------------------------------------- (b) leading-order state at t = 0
# Local frame tau_0 = (1,0), n_0 = (0,1) = perp(tau_0).  Ts = Q0[sigma].tau_0 (our T_0 sigma_0) and
# Qn = Q0[sigma].n_0 (= -U_0) are free; e2 = (tau2, n2) in this frame satisfies
#   n2 - A Qn = -sigma/2                     (proof of LM26 Prop. 7.4; our n_{0,2} + A U_0 = -sigma_0/2),
#   tau2 - A Ts = -(1-A^2) gamma/2           (tangential identity: our Prop. 2.10 = LM26 (3.7b), left relation).
sig, gam, kap, Ts, Qn, gn = sp.symbols('sigma gamma kappa T_s Q_n g_n', real=True)
c = (1 - A**2)*kap
tau2 = A*Ts - (1 - A**2)*gam/2
n2 = A*Qn - sig/2
E2 = sp.Matrix([tau2, n2])
Gp = sp.Matrix([Ts, Qn + A*(c + sig/2)])               # grad p at leading order, LM26 (7.11)
Gg = sp.Matrix([c + sig/4, gn])                         # grad g, LM26 (7.9): tangential, i.e. g_n = 0 for LM26
U = (-Gp + A*perp(Gg) - (A/2)*E2)/(1 + A**2/2)          # LM26 (4.10)
M_ = A*U - perp(Gg) - (A*U + E2)/2                      # LM26 (4.6)
V = U + A*M_
assert all(zero(x) for x in V - (-Gp - A*E2))
DL = 2*M_ + 2*E2                                        # D of (3.5) at rho = A
aL = E2 + (A/(1 + A**2))*V
W = A*U + E2
delta_star = -(1 - A**2)*gam/(2*(1 + A**2))            # optimal tilt, Lemma 5.2 with rho_m = A
delta_LM = -(1 - A**2)*gam/(2 + A**2)
assert zero(W[1] + sig/2)                               # (Au+e2).n_0 = -sigma/2, LM26 (7.10)
assert zero(W[0].subs(gn, 0) - delta_LM)                # (Au+e2).tau_0 = delta_LM
assert zero(DL[0].subs(gn, 0) - delta_LM)               # D.tau_0 = delta_LM
assert zero(delta_LM - 2*(1 + A**2)/(2 + A**2)*delta_star)
assert zero(DL[0] - (-(1 - A**2)*gam + 4*gn)/(2 + A**2))   # tilt as a function of the normal part g_n of grad g
assert zero(DL[0].subs(gn, (1 - A**2)*gam/4))           # delta = 0 needs g_n = (1-A^2)gamma/4 (not LM26's layer)
assert zero(DL[1] - (-sig - 2*c))                       # D.n_0 = -sigma - 2c = F_0/|d_alpha z_0| (rho_m = A)
# a agrees with Lemma 5.2 of the manuscript at rho_m = A
rm = sp.symbols('rho_m', real=True); mum = 1 + A*rm
a_n = -sig/2 - A**2*(1 - rm**2)*kap/mum
a_t = -(1 - A**2)*gam/(2*mum)
assert zero(aL[1] - a_n.subs(rm, A)) and zero(aL[0] - a_t.subs(rm, A))
# exact hull of LM26's unchanged layer at t = 0 (g_n = 0)
hullL = sp.expand(dot(DL, DL - 2*aL).subs(gn, 0))
assert zero(hullL - (2*c*(2*c + sig) - delta_LM**2)/(1 + A**2))
assert zero(hullL*(1 + A**2)/(2*(1 - A**2))
            - (kap*(2*(1 - A**2)*kap + sig) - (1 - A**2)*gam**2/(2*(2 + A**2)**2)))
# general identity of Lemma 5.2 (any rho_m, any tilt delta), with D_n, a_n, a_tau of (5.4)
delta = sp.symbols('delta', real=True)
D_n = -sig - 2*(1 - A*rm)*kap
Phi = kap*(2*(1 - A*rm)*kap + sig) - (1 - A**2)*gam**2/(8*mum)
ident = D_n*(D_n - 2*a_n) + delta*(delta - 2*a_t) - (2*(1 - A**2)/mum*Phi + (delta - a_t)**2)
assert zero(ident)
assert zero(hullL - (2*(1 - A**2)/mum*Phi + (delta - a_t)**2).subs({rm: A, delta: delta_LM}))
assert zero((delta_LM - delta_star) - A**2/(2 + A**2)*delta_star)     # extra term (A^2 delta_*/(2+A^2))^2
# N = 1 multilayer normalisation of LM26 (8.38), (8.43): D.tau_0 -> 2 delta_*, same margin as delta = 0
assert zero((delta*(delta - 2*a_t)).subs(delta, 2*a_t) - (delta*(delta - 2*a_t)).subs(delta, 0))
# gain of the optimal tilt on the gamma^2 coefficient (rho_m = A)
ratio = sp.simplify(((1 - A**2)/(8*(1 + A**2)))/((1 - A**2)/(2*(2 + A**2)**2)))   # ours / LM26's layer
assert zero(ratio - (2 + A**2)**2/(4*(1 + A**2))) and zero(ratio - (1 + A**4/(4*(1 + A**2))))
assert zero(sp.Rational(9, 8) - ratio - (1 - A**2)*(1 + 2*A**2)/(8*(1 + A**2)))   # ratio < 9/8 for |A| < 1
print("(b) LM26's single layer at t = 0, rho_m = A, c = (1-A^2) kappa, grad g tangential (LM26 (7.9)):")
print("    D.tau_0 = (Au+e2).tau_0 = -(1-A^2) gamma/(2+A^2) = delta_LM = [2(1+A^2)/(2+A^2)] delta_*;"
      "  D.n_0 = -sigma - 2c;  a = a^(0) of Lemma 5.2")
print("    strict hull at t = 0  <=>  kappa(2(1-A^2) kappa + sigma) < (1-A^2) gamma^2 / (2 (2+A^2)^2)")
print("    D.(D-2a) = (2(1-A^2)/mu_m) Phi(kappa) + (delta - delta_*)^2 (Lemma 5.2);"
      " delta_LM - delta_* = A^2 delta_*/(2+A^2)")
print("    gain of delta_* on the gamma^2 term: (2+A^2)^2/(4(1+A^2)) = 1 + A^4/(4(1+A^2)) in [1, 9/8)")

# ---------------------------------------------------------------- (c) numerics on the overturned example
# Periodic Nystrom discretisation of the adjoint Neumann-Poincare operator K, as in check_vorticity.py.
N = 512; h = 2*np.pi/N; al = np.arange(N)*h
z = (al - 1.5*np.sin(al)) + 1j*(0.8*np.sin(al)); dz = (1 - 1.5*np.cos(al)) + 1j*(0.8*np.cos(al))
d2z = (1.5*np.sin(al)) + 1j*(-0.8*np.sin(al))
L = np.abs(dz); tau = dz/L; nrm = 1j*tau
ZA, ZB = np.meshgrid(z, z, indexing='ij'); DA, DB = np.meshgrid(dz, dz, indexing='ij')
with np.errstate(all='ignore'):
    k = (np.abs(DB)/(4*np.pi*np.abs(DA)))*np.imag(DA/np.tan((ZB - ZA)/2))
np.fill_diagonal(k, -np.imag(d2z/dz)/(4*np.pi)); K = k*h; I = np.eye(N)
expected = {0.0: (0.2350, 0.2350), 0.3: (0.2297, 0.2294), -0.6: (0.2265, 0.2224), 0.9: (0.2959, 0.2778)}
print("(c) overturned example z_0 = (alpha - 1.5 sin alpha, 0.8 sin alpha), rho_m = A, minima over T:")
for Av, (kx, kl) in expected.items():
    s = np.linalg.solve(I - 2*Av*K, -2*nrm.imag)        # (I - 2AK) sigma = -2 n_2
    g = np.linalg.solve(I + 2*Av*K, -2*tau.imag)        # (I + 2AK) gamma = -2 tau_2
    B = 1 - Av**2; q = 1 - Av*Av; muv = 1 + Av*Av
    kmax = (-s + np.sqrt(s**2 + q*B*g**2/muv))/(4*q)                    # root of Phi
    dLM = -B*g/(2 + Av**2)
    klm = (-s + np.sqrt(s**2 + 4*dLM**2))/(4*B)                         # LM26's layer, exact hull
    kst = np.maximum(-s, 0)/(2*B)                                       # LM26's stated range cbar < 1/2
    st = s >= 0
    assert np.all(klm <= kmax + 1e-12)
    assert round(kmax.min(), 4) == kx and round(klm.min(), 4) == kl
    assert np.all(kst[st] == 0) and np.all(klm[st] > 0.2)
    print(f"    A={Av:+.1f}: min kappa_max = {kmax.min():.4f}, min kappa_LM (LM26 layer, exact hull) = {klm.min():.4f}, "
          f"min kappa_LM/kappa_max = {np.min(klm/kmax):.4f}; sigma_0 >= 0 on {st.mean():.1%} of the parameter "
          f"interval, {np.sum(L[st])/np.sum(L):.1%} of the arclength, min|gamma_0| there = {np.abs(g[st]).min():.3f}; "
          f"LM26's stated range is empty there")
    if Av == 0.9:
        assert np.min(klm/kmax) > 0.937
# scalar example: A = 0.5, sigma_0 = -1, gamma_0 = 1
Av, s, g = 0.5, -1.0, 1.0; B = 1 - Av**2
stated = -s/(2*B); exact = (-s + np.sqrt(s*s + 4*(B*g/(2 + Av**2))**2))/(4*B)
kx = (-s + np.sqrt(s*s + B*B*g*g/(1 + Av**2)))/(4*B)
assert (round(stated, 4), round(exact, 4), round(kx, 4)) == (0.6667, 0.7340, 0.7347)
print(f"    scalar example A=0.5, sigma_0=-1, gamma_0=1: LM26 stated {stated:.4f}, LM26 exact {exact:.4f}, "
      f"kappa_max {kx:.4f}")

# ---------------------------------------------------------------- (d) M versus d_alpha G_+ - d_alpha G_-
# Boundary fluxes F = D.(d_alpha z)^perp (the flux of D through the boundary per unit alpha, as in Sec. 4).
# (d1) pointwise, with D = w - 2 grad^perp g (part (a)) and LM26's d_alpha G = grad g . d_alpha z:
#      D.(d_alpha z)^perp = w.(d_alpha z)^perp - 2 d_alpha G.
T = sp.Matrix([p1, p2])
assert zero(dot(w - 2*gperp, perp(T)) - (dot(w, perp(T)) - 2*dot(gradg, T)))
# (d2) for any divergence-free w = grad^perp psi and the strip map X = Z + (lambda + A) nu, lambda in [-1,1]:
#      w(z+).(d_alpha z+)^perp - w(z-).(d_alpha z-)^perp = d_alpha int_{-1}^{1} w(X).nu^perp dlambda,
#      and nu = t kappa n_0 gives nu^perp = -t kappa tau_0.  Checked on explicit data.
al_, lam, t, x1, x2 = sp.symbols('alpha lambda t x1 x2', real=True)
psi = x1**2*x2 - x2**3/3 + 2*x1*x2 + x2                 # stream function of w (any polynomial)
gfun = x1*x2**2 + 3*x2 - x1**3/5                        # any potential g
def wf(P): return sp.Matrix([-sp.diff(psi, x2), sp.diff(psi, x1)]).subs({x1: P[0], x2: P[1]}, simultaneous=True)
def gradgf(P): return sp.Matrix([sp.diff(gfun, x1), sp.diff(gfun, x2)]).subs({x1: P[0], x2: P[1]}, simultaneous=True)
th0 = sp.sin(al_)/2                                     # tau_0 = (cos th0, sin th0), n_0 = perp(tau_0)
tau0 = sp.Matrix([sp.cos(th0), sp.sin(th0)]); n0 = perp(tau0)
kapf = 1 + sp.cos(al_)/3
Z = sp.Matrix([al_ + t*sp.sin(al_), sp.cos(al_) + t**2*sp.sin(2*al_)])   # any central curve
X = Z + (lam + A)*t*kapf*n0
zp, zm = X.subs(lam, 1), X.subs(lam, -1)
flux = lambda P: dot(wf(P), perp(sp.diff(P, al_)))
bracket = flux(zp) - flux(zm)
inner = sp.integrate(sp.expand(dot(wf(X), tau0)), (lam, -1, 1))
assert zero(bracket + sp.diff(t*kapf*inner, al_))
# undifferentiated form: psi(z+) - psi(z-) = -t kappa int_{-1}^{1} w(X).tau_0 dlambda, so that a stream function
# Psi = psi - 2g of D = w - 2 grad^perp g has Psi(z+) - Psi(z-) = -t kappa int w(X).tau_0 dlambda - 2(G_+ - G_-)
psi_at = lambda P: psi.subs({x1: P[0], x2: P[1]}, simultaneous=True)
assert zero(psi_at(zp) - psi_at(zm) + t*kapf*inner)
# D = w - 2 grad^perp g on the same strip:  M = F+ - F- = bracket - 2 (d_alpha G+ - d_alpha G-)
Dflux = lambda P: dot(wf(P) - 2*perp(gradgf(P)), perp(sp.diff(P, al_)))
Gp_, Gm_ = gfun.subs({x1: zp[0], x2: zp[1]}, simultaneous=True), gfun.subs({x1: zm[0], x2: zm[1]}, simultaneous=True)
assert zero(Dflux(zp) - Dflux(zm) - (bracket - 2*(sp.diff(Gp_, al_) - sp.diff(Gm_, al_))))
# (d3) at order t^0 the bracket vanishes; at order t it is -2 d_alpha(kappa w.tau_0) at the initial curve
ser = sp.expand(bracket)                                # a polynomial in t
assert zero(ser.coeff(t, 0))
assert zero(ser.coeff(t, 1) + 2*sp.diff(kapf*dot(wf(Z.subs(t, 0)), tau0), al_))
print("(d) with F^pm = D.(d_alpha z_pm)^perp and M = F^+ - F^- on LM26's layer (rho_m = A):")
print("    F^pm = (A U^pm + e2).(d_alpha z_pm)^perp - 2 d_alpha G_pm  (exact),")
print("    Psi^+ - Psi^- = -t kappa int_{-1}^{1} (Au+e2)(X).tau_0 dlambda - 2 (G_+ - G_-)  (exact; Psi stream function of D),")
print("    M = d_alpha int_{-1}^{1} (Au+e2)(X).(t kappa n_0)^perp dlambda - 2 (d_alpha G_+ - d_alpha G_-)  (exact),")
print("      = -2t d_alpha(kappa (Au+e2).tau_0) - 2 (d_alpha G_+ - d_alpha G_-) + O(t^2)  (smooth fields).")
print("    So M and d_alpha G_+ - d_alpha G_- agree (up to the factor -2) only at order t^0. LM26 cancel the O(t) part")
print("    of d_alpha G_+ - d_alpha G_- (Lemma 5.3); with (Au+e2).tau_0 = delta_LM + O(t) from (b) this gives")
print("    M = -2t d_alpha(kappa delta_LM) + O(t^2) on their layer, not O(t^2); and, since d_lambda G = (G_+ - G_-)/2")
print("    = O(t^2) (LM26 Prop. 7.2), Psi^+ - Psi^- = -2t kappa delta_LM + O(t^2) there.")
print("check_lm26_layer.py: all assertions passed")
