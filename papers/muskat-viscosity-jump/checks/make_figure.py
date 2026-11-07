import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
out_pdf, out_png = sys.argv[1], sys.argv[2]
N=512; h=2*np.pi/N; al=np.arange(N)*h
A=0.3; rm=A; mu=1+A*rm; t=0.25; theta=0.8
z=(al-1.5*np.sin(al))+1j*(0.8*np.sin(al)); dz=(1-1.5*np.cos(al))+1j*(0.8*np.cos(al)); d2z=(1.5*np.sin(al))+1j*(-0.8*np.sin(al))
L=np.abs(dz); tau=dz/L; n=1j*tau
ZA,ZB=np.meshgrid(z,z,indexing='ij'); DA,DB=np.meshgrid(dz,dz,indexing='ij')
with np.errstate(all='ignore'):
    k=(np.abs(DB)/(4*np.pi*np.abs(DA)))*np.imag(DA/np.tan((ZB-ZA)/2))
np.fill_diagonal(k,-np.imag(d2z/dz)/(4*np.pi)); K=k*h; I=np.eye(N)
sig=np.linalg.solve(I-2*A*K,-2*n.imag); gam=np.linalg.solve(I+2*A*K,-2*tau.imag); U0=-K@sig
kmax=(-sig+np.sqrt(sig**2+(1-A*rm)*(1-A**2)*gam**2/mu))/(4*(1-A*rm))
kold=np.maximum(-sig,0)/(2*(1-A*rm))
kap=theta*kmax; cp,cm=(1+rm)*kap,(1-rm)*kap
zc=z+t*U0*n; zp=zc+t*cp*n; zmn=zc-t*cm*n
R=N//2
roll=lambda x: np.roll(x,R)
al_pl=np.roll(al,R); al_pl=np.where(al_pl>=np.pi,al_pl-2*np.pi,al_pl)
z,zp,zmn,sig,gam,kmax,kold=[roll(x) for x in (z,zp,zmn,sig,gam,kmax,kold)]
z=np.where(al_pl<0,z-2*np.pi,z); zp=np.where(al_pl<0,zp-2*np.pi,zp); zmn=np.where(al_pl<0,zmn-2*np.pi,zmn)
al=al_pl
INK,MUTED,GRID='#0b0b0b','#52514e','#e6e5e0'; BLUE,ORANGE,AQUA='#2a78d6','#eb6834','#1baf7a'
plt.rcParams.update({'font.size':9.5,'axes.edgecolor':MUTED,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'font.family':'serif'})
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(7.6,4.2),gridspec_kw={'width_ratios':[1.25,1]})
# panel (a)
stable=sig>=0
def segs(mask):
    idx=np.where(mask)[0]; 
    if len(idx)==0: return []
    breaks=np.where(np.diff(idx)>1)[0]; parts=np.split(idx,breaks+1); return parts
ax1.fill(np.concatenate([zp.real,zmn.real[::-1]]),np.concatenate([zp.imag,zmn.imag[::-1]]),color=AQUA,alpha=0.25,lw=0,label=r'mixing zone at $t=0.25$ (first order in $t$)')
ax1.plot(zp.real,zp.imag,color=MUTED,lw=0.8); ax1.plot(zmn.real,zmn.imag,color=MUTED,lw=0.8)
first=True
for p in segs(~stable):
    ax1.plot(z.real[p],z.imag[p],color=BLUE,lw=2,label=r'$\sigma_0<0$ (Saffman–Taylor unstable)' if first else None); first=False
first=True
for p in segs(stable):
    ax1.plot(z.real[p],z.imag[p],color=ORANGE,lw=2,label=r'$\sigma_0\geq0$, $\gamma_0\neq0$ (stable, vortex sheet)' if first else None); first=False
ax1.set_aspect('equal'); ax1.set_xlabel(r'$x_1$'); ax1.set_ylabel(r'$x_2$')
ax1.grid(color=GRID,lw=0.6); ax1.set_axisbelow(True)
ax1.legend(loc='upper center',bbox_to_anchor=(0.5,-0.3),fontsize=7.5,frameon=False,ncol=1)
ax1.set_title(r'(a) overturned interface, $A=0.3$, $\rho_{\rm m}=A$',fontsize=9,color=INK,loc='left')
# panel (b)
for p in segs(stable):
    ax2.axvspan(al[p[0]],al[p[-1]],color=GRID,lw=0)
ax2.plot(al,kmax,color=BLUE,lw=2,label=r'$\kappa_{\max}$ (this paper)')
ax2.plot(al,kold,color=ORANGE,lw=2,label=r'$\max(-\sigma_0,0)/(2(1-A\rho_{\rm m}))$')
ax2.set_xlim(-np.pi,np.pi); ax2.set_ylim(bottom=0); ax2.set_xlabel(r'$\alpha$'); ax2.set_ylabel('admissible opening function')
ax2.text(al[np.where(stable)[0][len(np.where(stable)[0])//2]],ax2.get_ylim()[1]*0.92,r'$\sigma_0\geq0$',ha='center',color=MUTED,fontsize=8)
ax2.grid(color=GRID,lw=0.6,axis='y'); ax2.set_axisbelow(True)
ax2.legend(loc='upper center',bbox_to_anchor=(0.5,-0.3),fontsize=7.5,frameon=False,ncol=1)
ax2.set_title(r'(b) upper bounds for $\kappa(\alpha)$',fontsize=9,color=INK,loc='left')
for ax in (ax1,ax2):
    for s in ('top','right'): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(out_pdf); fig.savefig(out_png,dpi=110)
print("stable fraction %.2f, min kmax %.3f, max kmax %.3f, min|gamma| on stable %.3f"%(stable.mean(),kmax.min(),kmax.max(),np.abs(gam[stable]).min()))
