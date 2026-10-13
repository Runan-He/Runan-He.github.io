#!/usr/bin/env python3
"""Figure for Section 5: spectral gap of the ODE blowup profile in the unified
power family  u_tt - Lap u = sgn(p-1)|u|^{p-1}u,  s = 2/(p-1).

  omega0(p)  = min(1, 2s+2)        true gap of L (eigenvalue-determined)
  2s+2       = -lambda_*           top of the second eigenvalue family (l=0, n=0)
  min(1, s)  = Ostermann's rate    (decay rate in Ostermann 2024, p>1; the free
                                    semigroup has the eigenvalue -s from u = const)
Markers: gaps found numerically by the monodromy test (power_spectrum.py), d=3.

Output: ../figures/power_gap.pdf
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'figures', 'power_gap.pdf')

# reference categorical palette (light mode), slots 1-3
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#d9d8d4'

plt.rcParams.update({
    'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
    'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.linewidth': 0.6,
    'font.family': 'serif', 'mathtext.fontset': 'cm',
})

s_of = lambda p: 2.0/(p - 1.0)
gap = lambda p: np.minimum(1.0, 2*s_of(p) + 2)

# numerically confirmed gaps (monodromy test, d = 3, see power_spectrum.py)
num = {3.0: 1.0, 7.0: 1.0, -2.0: 2/3, -3.0: 1.0}

fig, (axL, axR) = plt.subplots(1, 2, figsize=(6.4, 2.7), sharey=True,
                               gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.08})

# ---- MEMS range p = -q < -1
pL = np.linspace(-12, -1.0005, 1500)
axL.plot(pL, 2*s_of(pL) + 2, color=MUTED, lw=1.0, ls='--', label=r'$2s+2$ (2nd family, $\ell=n=0$)')
axL.plot(pL, gap(pL), color=BLUE, lw=2.0, label=r'$\omega_0=\min\{1,2s+2\}$')
axL.axvline(-3, color=GRID, lw=0.8, zorder=0)
axL.text(-3.15, 0.08, r'$q=3$', ha='right', color=MUTED, fontsize=8)
axL.set_xlim(-12, -1)
axL.set_xlabel(r'$p=-q$  (MEMS range, $s=-\sigma\in(-1,0)$)')
axL.set_ylabel(r'decay exponent')

# ---- power range p > 1
pR = np.linspace(1.0005, 12, 1500)
axR.plot(pR, gap(pR), color=BLUE, lw=2.0)
axR.plot(pR, np.minimum(1.0, s_of(pR)), color=ORANGE, lw=1.6, ls=':',
         label=r'$\min\{1,s\}$ (rate in [Ostermann 2024])')
axR.plot(pR, 2*s_of(pR) + 2, color=MUTED, lw=1.0, ls='--')
axR.axvline(3, color=GRID, lw=0.8, zorder=0)
axR.text(3.2, 0.08, r'$p=3$', ha='left', color=MUTED, fontsize=8)
axR.set_xlim(1, 12)
axR.set_xlabel(r'$p>1$  (power range, $s=2/(p-1)>0$)')

for ax in (axL, axR):
    ax.set_ylim(0, 2.5)
    ax.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for sp_ in ('top', 'right'):
        ax.spines[sp_].set_visible(False)
for p, g in num.items():
    ax = axL if p < 0 else axR
    ax.plot([p], [g], 'o', ms=6, mfc=AQUA, mec='white', mew=1.2, zorder=5,
            label='numerical gap ($d=3$)' if p == -2.0 else None)

# one legend for the figure (all series named once)
h1, l1 = axL.get_legend_handles_labels()
h2, l2 = axR.get_legend_handles_labels()
fig.legend(h1 + h2, l1 + l2, loc='upper center', ncol=2, frameon=False, fontsize=8,
           bbox_to_anchor=(0.5, 1.06))
fig.savefig(OUT, bbox_inches='tight')
print('wrote', os.path.normpath(OUT))
