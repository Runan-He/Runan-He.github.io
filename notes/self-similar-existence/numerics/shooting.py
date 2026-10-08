"""Two-point shooting for the smooth global self-similar profiles (Section 6 of the note).
Usage: python3 shooting.py <n> [eta_m]   (run inner_phase.py <n> first, it writes phases_n<n>.npy)
"""
import numpy as np, sys, warnings
from profile_ode import *
warnings.filterwarnings('ignore')
eta_m = float(sys.argv[2]) if len(sys.argv) > 2 else 0.85
def fwd_pt(c):
    s = forward(c, eta_m); v, q = s.y[:, -1]; return v, q/eta_m          # (v, v') at eta_m
def bwd_pt(c1):
    s = backward(c1, eta_m); v, q = s.y[:, -1]; return v, q/eta_m
# monotonicity of v_b(eta_m; c1) in c1 on a bracket around a
c1grid = np.linspace(0.3, 3.0, 109); vb = np.array([bwd_pt(x)[0] for x in c1grid])
ia = int(np.argmin(abs(c1grid-a))); lo = ia; hi = ia
while lo > 0 and vb[lo-1] < vb[lo]: lo -= 1
while hi < len(c1grid)-1 and vb[hi+1] > vb[hi]: hi += 1
BR = (c1grid[lo], c1grid[hi])
print(f"n={n}, eta_m={eta_m}: v_b(eta_m;c1) increasing on c1 in [{BR[0]:.3f},{BR[1]:.3f}] with values [{vb[lo]:.3f},{vb[hi]:.3f}]; a={a:.4f}")
def g(c):
    vf, pf = fwd_pt(c)
    c1 = brentq(lambda x: bwd_pt(x)[0]-vf, BR[0], BR[1], xtol=1e-14)
    return pf - bwd_pt(c1)[1], c1
cs = np.logspace(-3.3, 0.25, 420); gs = []
for c in cs:
    try: gs.append(g(c)[0])
    except ValueError: gs.append(np.nan)
gs = np.array(gs); print(f"  scan: {np.sum(np.isnan(gs))} of {len(cs)} values of c could not be matched (v_f outside the monotone range)")
roots = []
for i in range(len(cs)-1):
    if np.isfinite(gs[i]) and np.isfinite(gs[i+1]) and gs[i]*gs[i+1] < 0:
        r = brentq(lambda c: g(c)[0], cs[i], cs[i+1], xtol=1e-14); roots.append(r)
roots = sorted(roots, reverse=True)
ph = np.load(f'phases_n{n}.npy'); A1, A2, B1, B2 = ph
print(f"  predicted ratio exp(-2pi/(3b)) = {np.exp(-2*np.pi/(3*b)):.5f}")
out = []
for k, c in enumerate(roots, 1):
    c1 = g(c)[1]
    sb = beyond_one(c1, 1e4); vend = sb.y[0, -1]; eend = np.exp(sb.t[-1])
    sf = forward(c, 0.999, dense=True); et = np.exp(sf.t); vv = sf.y[0]
    ncross = np.sum(np.diff(np.sign(vv - a*et**(2/3))) != 0)
    mono = np.all(np.diff(vv) >= -1e-12) and np.all(np.diff(sb.y[0]) >= -1e-12)
    ratio = roots[k]/c if k < len(roots) else float('nan')
    cpred = np.exp((2/(3*b))*(A2 - B2 - np.pi*k))   # formal matching with pi k (phase mod pi)
    print(f"  k={k}: c_k={c:.8f}  v_k(1)={c1:.6f}  ratio c_(k+1)/c_k={ratio:.5f}  b_k=lim v/eta^(2/3)~{vend/eend**(2/3):.5f} (eta={eend:.0f}) crossings of a*eta^(2/3) on (0,1): {ncross}  monotone: {mono}  c_k(formal,pi k)={cpred:.6f}")
    out.append([k, c, c1, vend/eend**(2/3), ncross])
np.save(f'roots_n{n}.npy', np.array(out))
