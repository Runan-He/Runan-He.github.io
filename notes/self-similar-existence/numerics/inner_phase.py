"""Inner phase (A1, A2) of V_0 and phase (B1, B2) of the regular linear mode at the cone.
Usage: python3 inner_phase.py <n>
"""
import numpy as np, sys
from profile_ode import *
def rhs_in(s, y):
    V, q = y; e2 = np.exp(2*s); return [q, q - (n-1)*q + e2/V**2]
r0 = 1e-4; V2 = 1/(2*n)
sol = solve_ivp(rhs_in, [np.log(r0), np.log(2e3)], [1+V2*r0**2, 2*V2*r0**2], method='DOP853', rtol=1e-13, atol=1e-15, dense_output=True)
lo, hi = (30, 2e3) if n <= 3 else (8, 300)
ss_ = np.linspace(np.log(lo), np.log(hi), 4000)
Z = sol.sol(ss_)[0]*np.exp(-2/3*ss_) - a
w = np.exp(-(1/3-n/2)*ss_)*Z
# include the next-order correction r^{2(1/3-n/2)} (cos 2.., sin 2.., const) as nuisance terms
M = np.vstack([np.cos(b*ss_), np.sin(b*ss_), np.exp((1/3-n/2)*ss_)*np.cos(2*b*ss_), np.exp((1/3-n/2)*ss_)*np.sin(2*b*ss_), np.exp((1/3-n/2)*ss_)]).T
coef, res, _, _ = np.linalg.lstsq(M, w, rcond=None)
A1 = np.hypot(coef[0], coef[1]); A2 = np.arctan2(-coef[1], coef[0])
print(f"n={n}: A1={A1:.6f} A2={A2:.6f} relres={np.sqrt(res[0])/np.linalg.norm(w):.1e}; V0 increasing: {np.all(np.diff(sol.y[0])>0)}")


# Linearised equation about a*eta^(2/3): the solution regular at eta=1, in s=ln(eta):
# (1-e^{2s}) W_ss + ((n-2) + e^{2s}/3) W_s + ((4/3)(n-4/3) + (2/9)e^{2s}) W = 0, W(1)=1
def rhs_lin(s, y):
    W, q = y; e2 = np.exp(2*s)
    return [q, -(((n-2)+e2/3)*q + (4/3*(n-4/3) + 2/9*e2)*W)/(1-e2)]
# Taylor series at eta=1 for the regular solution: at x=eta-1, leading order gives W'(1) = -(4/3)(n-4/3)+2/9 ... use the ODE
# Simpler: start at eta=1-d with a few Taylor terms computed numerically by series recursion.
import numpy.polynomial.polynomial as P
Mt = 20; t = np.zeros(Mt+1); t[0] = 1.0
A = np.zeros(Mt+2); A[1] = -2; A[2] = -1
B = np.array([(n-1)*(-1)**m for m in range(Mt+2)], float); B[0] -= 2/3; B[1] -= 2/3
Cc = np.array([2/9] + [4/3*(n-4/3)*(m+1)*(-1)**m for m in range(0, Mt+1)], float)  # 2/9 + (4/3)(n-4/3)/(1+x)^2
Cc = np.array([2/9 + 4/3*(n-4/3)] + [4/3*(n-4/3)*(m+1)*(-1)**m for m in range(1, Mt+2)], float)
for m in range(Mt):
    def resid(tm1):
        tt = t.copy(); tt[m+1] = tm1; v = tt[:m+2]
        tot = P.polyadd(P.polyadd(P.polymul(A, P.polyder(v,2)), P.polymul(B, P.polyder(v))), P.polymul(Cc, v))
        return tot[m]
    r0_, r1_ = resid(0.0), resid(1.0); t[m+1] = -r0_/(r1_-r0_)
d = 0.02; x = -d
W = sum(t[m]*x**m for m in range(Mt+1)); Wp = sum(m*t[m]*x**(m-1) for m in range(1, Mt+1))
sl = solve_ivp(rhs_lin, [np.log(1-d), np.log(1e-7)], [W, (1-d)*Wp], method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
ss2 = np.linspace(np.log(1e-7), np.log(1e-4), 3000)
w2 = sl.sol(ss2)[0]*np.exp(-(1-n/2)*ss2)
M2 = np.vstack([np.cos(b*ss2), np.sin(b*ss2)]).T
coef2, res2, _, _ = np.linalg.lstsq(M2, w2, rcond=None)
B1 = np.hypot(*coef2); B2 = np.arctan2(-coef2[1], coef2[0])
print(f"n={n}: regular-at-1 linear mode: amplitude B1={B1:.6f}, phase B2={B2:.6f}, rel. fit residual={np.sqrt(res2[0])/np.linalg.norm(w2):.2e}")
print(f"theta = B2 - A2 = {B2-A2:.6f}")
print(f"predicted spacing of ln c_k: pi*2/(3b) = {2*np.pi/(3*b):.6f}  -> ratio c_(k+1)/c_k = {np.exp(-2*np.pi/(3*b)):.6f}")


np.save(f'phases_n{n}.npy', np.array([A1, A2, B1, B2]))
