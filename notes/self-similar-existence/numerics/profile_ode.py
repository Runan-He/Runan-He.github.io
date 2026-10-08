"""Numerical sanity checks for the self-similar ODE
   (1-eta^2) v'' + ((n-1)/eta - (2/3) eta) v' + (2/9) v - 1/v^2 = 0 ,  v(0)=c, v'(0)=0.
Checks: (i) inner profile V0 and its log-oscillation around a*r^{2/3};
        (ii) the regular-at-1 solution of the linearised equation and its phase B2 as eta->0;
        (iii) two-point shooting for the smooth global solutions v_k, and the ratio c_{k+1}/c_k.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve
import sys

n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
cstar = (9/2)**(1/3)
a = (2/3*(n-4/3))**(-1/3)
Delta = (n/2-7/3)**2 - 8/3
b = np.sqrt(-Delta) if Delta < 0 else float('nan')
alpha = n/2 + 1/6            # singular exponent at the light cone eta=1

def rhs(eta, y):
    v, p = y
    return [p, (-((n-1)/eta - 2/3*eta)*p - 2/9*v + 1/v**2)/(1-eta**2)]

def rhs_s(s, y):               # same ODE in s = ln(eta): y = (v, dv/ds)
    v, q = y; e2 = np.exp(2*s)
    # (1-e2)(q_s - q) + ((n-1) - 2/3 e2) q + (2/9 e2) v - e2/v^2 = 0
    qs = q - (((n-1) - 2/3*e2)*q + 2/9*e2*v - e2/v**2)/(1-e2)
    return [q, qs]

def forward(c, eta_end, dense=False):
    """regular solution with v(0)=c, integrated in s=ln(eta) from a tiny eta."""
    v2 = (1/c**2 - 2/9*c)/(2*n)
    eta0 = 1e-4*min(c**1.5, 1.0)
    y0 = [c + v2*eta0**2, 2*v2*eta0**2]
    sol = solve_ivp(rhs_s, [np.log(eta0), np.log(eta_end)], y0, method='DOP853',
                    rtol=1e-12, atol=1e-14, dense_output=dense)
    return sol

def taylor_at_one(c1, M=24):
    """Taylor coefficients of the smooth solution at eta=1: v = sum_m t_m (eta-1)^m, v(1)=c1.
    Obtained by solving the ODE order by order (non-resonant: alpha not an integer)."""
    # write x = eta-1: (1-eta^2) = -x(2+x); (n-1)/eta - 2/3 eta expanded in x.
    # Use a numerical recursion via truncated power series arithmetic.
    import numpy.polynomial.polynomial as P
    t = np.zeros(M+1); t[0] = c1
    # coefficient series: A(x) = -x(2+x);  B(x) = (n-1)/(1+x) - 2/3 (1+x);  C = 2/9
    A = np.zeros(M+2); A[1] = -2; A[2] = -1
    B = np.array([(n-1)*(-1)**m for m in range(M+2)], float); B[0] -= 2/3; B[1] -= 2/3
    for m in range(0, M):
        # collect the coefficient of x^m in A v'' + B v' + C v - 1/v^2 using known t_0..t_m, unknown t_{m+1}
        # coefficient of x^m: from A v'': sum_j A[j]*(m-j+2)(m-j+1) t_{m-j+2}, j>=1 -> involves t_{m+1} (j=1) and lower
        #                      from B v'': sum_j B[j]*(m-j+1) t_{m-j+1} -> t_{m+1} (j=0) and lower
        def resid(tm1):
            tt = t.copy(); tt[m+1] = tm1
            vser = tt[:m+2]
            # 1/v^2 series up to x^m
            inv = np.zeros(m+1); inv[0] = 1/vser[0]
            for k in range(1, m+1):
                inv[k] = -sum(vser[j]*inv[k-j] for j in range(1, k+1))/vser[0]
            inv2 = P.polymul(inv, inv)[:m+1]
            d1 = P.polyder(vser); d2 = P.polyder(vser, 2)
            tot = P.polyadd(P.polyadd(P.polymul(A, d2), P.polymul(B, d1)), 2/9*vser)
            tot = P.polysub(tot, inv2)
            return tot[m] if m < len(tot) else 0.0
        r0 = resid(0.0); r1 = resid(1.0)
        t[m+1] = -r0/(r1-r0)
    return t

def backward(c1, eta_end, delta=0.02, dense=False):
    """smooth solution through eta=1 with v(1)=c1, started at eta=1-delta from its Taylor series,
    integrated (in s=ln eta) down to eta_end<1."""
    t = taylor_at_one(c1); x = -delta
    v = sum(t[m]*x**m for m in range(len(t))); vp = sum(m*t[m]*x**(m-1) for m in range(1, len(t)))
    eta1 = 1-delta
    sol = solve_ivp(rhs_s, [np.log(eta1), np.log(eta_end)], [v, eta1*vp], method='DOP853',
                    rtol=1e-12, atol=1e-14, dense_output=dense)
    return sol

def beyond_one(c1, eta_end, delta=0.02):
    t = taylor_at_one(c1); x = delta
    v = sum(t[m]*x**m for m in range(len(t))); vp = sum(m*t[m]*x**(m-1) for m in range(1, len(t)))
    eta1 = 1+delta
    def ev(eta, y): return y[0]-1e-6
    ev.terminal = True
    return solve_ivp(rhs_s, [np.log(eta1), np.log(eta_end)], [v, eta1*vp], method='DOP853',
                     rtol=1e-11, atol=1e-13, events=ev, dense_output=True)

# --- robust versions (small v(1): the Taylor series at the cone has a small radius of convergence)
def _taylor_data(c1, side):
    t = taylor_at_one(c1, M=24)
    good = np.isfinite(t); M = int(np.argmin(good)) - 1 if not good.all() else len(t) - 1
    t = t[:M+1]
    delta = 0.02
    while abs(t[M])*delta**M > 1e-13 and delta > 1e-4: delta /= 2
    x = side*delta
    v = sum(t[m]*x**m for m in range(M+1)); vp = sum(m*t[m]*x**(m-1) for m in range(1, M+1))
    return 1+x, v, vp

def backward(c1, eta_end, dense=False):
    eta1, v, vp = _taylor_data(c1, -1)
    return solve_ivp(rhs_s, [np.log(eta1), np.log(eta_end)], [v, eta1*vp], method='DOP853',
                     rtol=1e-12, atol=1e-14, dense_output=dense)

def beyond_one(c1, eta_end):
    eta1, v, vp = _taylor_data(c1, +1)
    def ev(eta, y): return y[0]-1e-6
    ev.terminal = True
    return solve_ivp(rhs_s, [np.log(eta1), np.log(eta_end)], [v, eta1*vp], method='DOP853',
                     rtol=1e-11, atol=1e-13, events=ev, dense_output=True)
