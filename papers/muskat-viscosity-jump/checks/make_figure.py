# Figure 1 (fig-turned.pdf) of the manuscript.
# Usage, from the paper folder:  python3 checks/make_figure.py fig-turned.pdf preview.png
#
# Overturned interface z0(a) = (a - 1.5 sin a, 0.8 sin a), A = 0.3, rho_m = A.
# sigma_0 and gamma_0 are computed by the periodic Nystrom (trapezoidal) discretisation of eq. (data), as in
# check_vorticity.py.  Panel (a): the initial interface z0 (its stable part {sigma_0 >= 0} highlighted) and the
# first-order strip z0 + t(U0 +- (1 +- rho_m) kappa) n0, kappa = 0.8 kappa_max, with its central curve z0 + t U0 n0,
# at the illustrative time t = 0.25.  Panel (b): three single-strip upper bounds for kappa along the interface:
#   kappa_max  -- the sharp range of the manuscript (optimal tilt delta_*),
#   kappa_LM   -- the single layer of LM26 with its exact hull condition (4.14) at t = 0, i.e. the range (x:LM-range):
#                 kappa_LM = (-sigma_0 + sqrt(sigma_0^2 + 4 T^2)) / (4 (1 - A^2)),  T = -(1 - A^2) gamma_0 / (2 + A^2),
#   the opening bound stated in LM26 (cbar < 1/2):  max(-sigma_0, 0) / (2 (1 - A^2)).
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
out_pdf, out_png = sys.argv[1], sys.argv[2]

N = 512; h = 2*np.pi/N; al = np.arange(N)*h
A = 0.3; rm = A; mu = 1 + A*rm; B = 1 - A**2; t = 0.25; theta = 0.8
z = (al - 1.5*np.sin(al)) + 1j*(0.8*np.sin(al)); dz = (1 - 1.5*np.cos(al)) + 1j*(0.8*np.cos(al))
d2z = (1.5*np.sin(al)) + 1j*(-0.8*np.sin(al))
L = np.abs(dz); tau = dz/L; n = 1j*tau
ZA, ZB = np.meshgrid(z, z, indexing='ij'); DA, DB = np.meshgrid(dz, dz, indexing='ij')
with np.errstate(all='ignore'):
    k = (np.abs(DB)/(4*np.pi*np.abs(DA)))*np.imag(DA/np.tan((ZB - ZA)/2))
np.fill_diagonal(k, -np.imag(d2z/dz)/(4*np.pi)); K = k*h; I = np.eye(N)
sig = np.linalg.solve(I - 2*A*K, -2*n.imag); gam = np.linalg.solve(I + 2*A*K, -2*tau.imag); U0 = -K@sig

kmax = (-sig + np.sqrt(sig**2 + (1 - A*rm)*B*gam**2/mu))/(4*(1 - A*rm))      # this paper, optimal tilt
T = -B*gam/(2 + A**2)                                                        # delta_LM = (Au+e2).tau0 on LM26's layer
klm = (-sig + np.sqrt(sig**2 + 4*T**2))/(4*B)                                # LM26 single layer, exact hull
kst = np.maximum(-sig, 0)/(2*B)                                              # LM26 stated (cbar < 1/2)

kap = theta*kmax; cp, cm = (1 + rm)*kap, (1 - rm)*kap
zc = z + t*U0*n; zp = zc + t*cp*n; zmn = zc - t*cm*n
stable = sig >= 0
# z0(alpha) lies outside the drawn strip when 0 is not between the normal offsets t(U0 - c-) and t(U0 + c+)
outside = (t*(U0 - cm) > 0) | (t*(U0 + cp) < 0)
report = ("A=%.1f, N=%d: sigma_0 >= 0 on %.4f of the parameter interval and %.4f of the arclength; "
          "min |gamma_0| there %.4f\n" % (A, N, stable.mean(), np.sum(L[stable])/np.sum(L), np.abs(gam[stable]).min())
          + "min kappa_max %.4f, min kappa_LM %.4f, min kappa_LM/kappa_max %.4f, stated bound on the stable set %.1e\n"
          % (kmax.min(), klm.min(), np.min(klm/kmax), kst[stable].max())
          + "at t=%.2f, z0 lies outside the drawn strip at %.0f%% of the nodes (%.0f%% of the stable nodes)"
          % (t, 100*outside.mean(), 100*outside[stable].mean()))

# plot over alpha in [-pi, pi)
R = N//2
roll = lambda x: np.roll(x, R)
al_pl = np.roll(al, R); al_pl = np.where(al_pl >= np.pi, al_pl - 2*np.pi, al_pl)
z, zc, zp, zmn, sig, kmax, klm, kst, stable = [roll(x) for x in (z, zc, zp, zmn, sig, kmax, klm, kst, stable)]
shift = np.where(al_pl < 0, -2*np.pi, 0.0)
z, zc, zp, zmn = z + shift, zc + shift, zp + shift, zmn + shift
al = al_pl

INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e6e5e0'
BLUE, ORANGE, GREEN = '#2a78d6', '#eb6834', '#1baf7a'
SHADE = '#e3e1da'          # marks {sigma_0 >= 0} in panel (b)
SHADE_A = '#a8a69e'        # marks {sigma_0 >= 0} on z0 in panel (a)
plt.rcParams.update({'font.size': 9.5, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'font.family': 'serif', 'mathtext.fontset': 'dejavuserif'})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.75), gridspec_kw={'width_ratios': [1.3, 1]})

def segs(mask):
    idx = np.where(mask)[0]
    if len(idx) == 0: return []
    return np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)

# ---------------- panel (a)
ax1.fill(np.concatenate([zp.real, zmn.real[::-1]]), np.concatenate([zp.imag, zmn.imag[::-1]]),
         color=BLUE, alpha=0.17, lw=0, zorder=2)
ax1.plot(zp.real, zp.imag, color=BLUE, lw=0.6, alpha=0.7, zorder=2)
ax1.plot(zmn.real, zmn.imag, color=BLUE, lw=0.6, alpha=0.7, zorder=2)
for p in segs(stable):
    ax1.plot(z.real[p], z.imag[p], color=SHADE_A, lw=6, solid_capstyle='butt', zorder=3)
ax1.plot(zc.real, zc.imag, color=BLUE, lw=1.3, zorder=3.5)
ax1.plot(z.real, z.imag, color=INK, lw=1.3, zorder=4)
ax1.set_aspect('equal'); ax1.set_xlabel(r'$x_1$'); ax1.set_ylabel(r'$x_2$')
ax1.set_xlim(-np.pi - 0.35, np.pi + 0.35); ax1.set_ylim(-1.35, 1.45)
ax1.grid(color=GRID, lw=0.6); ax1.set_axisbelow(True)
ax1.set_title(r'(a) interface and strip, $A=0.3$, $\rho_{\rm m}=A$', fontsize=9.5, color=INK, loc='left')
hand_a = [Line2D([], [], color=INK, lw=1.3, label=r'initial interface $z_0$'),
          Line2D([], [], color=SHADE_A, lw=6, solid_capstyle='butt', label=r'part of $z_0$ where $\sigma_0\geq0$'),
          Line2D([], [], color=BLUE, lw=1.3, label=r'central curve $z_0+tU_0n_0$'),
          Patch(facecolor=BLUE, alpha=0.17, edgecolor=BLUE, lw=0.6,
                label=r'strip for $\kappa=0.8\,\kappa_{\max}$ (first order in $t$)')]
ax1.legend(handles=hand_a, loc='upper left', bbox_to_anchor=(0.0, -0.31), fontsize=8.2, frameon=False, ncol=1,
           handlelength=2.4)
ax1.text(0.99, 0.03, r'$t=0.25$', transform=ax1.transAxes, ha='right', va='bottom', color=MUTED, fontsize=8)

# ---------------- panel (b)
for p in segs(stable):
    ax2.axvspan(al[p[0]], al[p[-1]], color=SHADE, lw=0, zorder=0)
ax2.plot(al, kst, color=ORANGE, lw=1.8, zorder=2)
ax2.plot(al, kmax, color=BLUE, lw=1.8, zorder=3)
ax2.plot(al, klm, color=GREEN, lw=1.8, dashes=(3.2, 2.4), zorder=4)
ax2.set_xlim(-np.pi, np.pi); ax2.set_ylim(0, 1.45)
ax2.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi]); ax2.set_xticklabels([r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'])
ax2.set_xlabel(r'$\alpha$'); ax2.set_ylabel(r'upper bound for $\kappa(\alpha)$')
st_idx = np.where(stable)[0]
ax2.text(al[st_idx[len(st_idx)//2]], 1.37, r'$\sigma_0\geq0$', ha='center', va='center', color=MUTED, fontsize=8)
ax2.grid(color=GRID, lw=0.6, axis='y'); ax2.set_axisbelow(True)
ax2.set_title(r'(b) single-strip bounds, $\rho_{\rm m}=A$', fontsize=9.5, color=INK, loc='left')
hand_b = [Line2D([], [], color=BLUE, lw=1.8, label=r'$\kappa_{\max}$ (optimal tilt $\delta_*$)'),
          Line2D([], [], color=GREEN, lw=1.8, dashes=(3.2, 2.4), label=r'single layer of [LM26], exact hull (4.14)'),
          Line2D([], [], color=ORANGE, lw=1.8, label=r'bound stated in [LM26] ($\bar c<\frac{1}{2}$)')]
ax2.legend(handles=hand_b, loc='upper left', bbox_to_anchor=(0.0, -0.24), fontsize=8.2, frameon=False, ncol=1,
           handlelength=2.4)
for ax in (ax1, ax2):
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
fig.tight_layout(w_pad=1.5)
fig.savefig(out_pdf, bbox_inches='tight', pad_inches=0.03); fig.savefig(out_png, dpi=150, bbox_inches='tight', pad_inches=0.03)

print(report)
