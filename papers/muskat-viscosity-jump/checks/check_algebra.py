import sympy as sp, numpy as np
A,rm,k,U,n2,P,Wp,Wm,vS,Vp,Vm,Q,Ap,Am_=sp.symbols('A rho_m kappa U n2 P Wp Wm vS Vp Vm Q Ap Am',real=True)
mu=1+A*rm
# Lemma 5.1: M0/|dz0| = 4(P-Wbar)/(1-rm^2)
Wbar=((1-rm)*Wp+(1+rm)*Wm)/2
vS0=P+A/2*((1+rm)*Wm-(1-rm)*Wp)
Fp=2/(1+rm)*(vS0-(1-A)*Wp)+2*n2; Fm=2/(1-rm)*((1+A)*Wm-vS0)+2*n2
print("M0 formula:", sp.simplify(Fp-Fm-4*(P-Wbar)/(1-rm**2))==0)
# ansatz: Wp=U+(1+rm)k, Wm=U-(1-rm)k, P=U ; sigma0=-2(n2+AU)
sub={Wp:U+(1+rm)*k, Wm:U-(1-rm)*k, P:U}
sig0=-2*(n2+A*U)
print("Wbar=U:", sp.simplify(Wbar.subs(sub)-U)==0)
print("vS0.n formula:", sp.simplify(vS0.subs(sub)-(mu*U-A*(1-rm**2)*k))==0)
print("F0 formula (+):", sp.simplify(Fp.subs(sub)-(-sig0-2*(1-A*rm)*k))==0, " (-):", sp.simplify(Fm.subs(sub)-(-sig0-2*(1-A*rm)*k))==0)
# Lemma 5.2
an=n2+A/mu*(mu*U-A*(1-rm**2)*k); Dn=-sig0-2*(1-A*rm)*k
print("a.n formula:", sp.simplify(an-(-sig0/2-A**2*(1-rm**2)*k/mu))==0)
print("hull at t=0:", sp.simplify(Dn*(Dn-2*an)-(-2*k*(1-A**2)/mu*Dn))==0)
# Lemma 4.2: RH on Gamma+ and Gamma-
dnp=2/(1+rm)*(vS-(1-A)*Vp); mS=((rm-A)*vS+sp.Rational(1,2)*(1-rm**2)*dnp)/(1-A**2)
mab=(vS+A*(1-rm)*Vp)/(1+A)
print("RH rho Gamma+:", sp.simplify(mab-mS-(1-rm)*Vp)==0, " u.n Gamma+:", sp.simplify((vS-A*mS)-mab)==0)
dnm=2/(1-rm)*((1+A)*Vm-vS); mS2=((rm-A)*vS+sp.Rational(1,2)*(1-rm**2)*dnm)/(1-A**2)
vb=vS-A*(1+rm)*Vm; mbe=-vb/(1-A)
print("RH rho Gamma-:", sp.simplify(mS2-mbe-(1+rm)*Vm)==0, " u.n Gamma-:", sp.simplify((vS-A*mS2)-vb/(1-A))==0)
# Lemma 4.3 mass identity
intFp=2/(1+rm)*(A*Q/2-(1-A*rm)*Ap); intFm=2/(1-rm)*((1-A*rm)*Am_+A*Q/2)
print("mass identity:", sp.simplify((intFp-intFm).subs(Q,(1-rm)*Ap+(1+rm)*Am_)+2*((1-rm)*Ap+(1+rm)*Am_)/(1-rm**2))==0)
# Example 3.6 flat
th,r=sp.symbols('theta r',real=True)
F=th*(1-r**2)/(1-A*r); m2=-F; v2=A*m2; mu_r=1+A*r
d2=2*((1-A**2)*m2-(r-A)*v2)/(1-r**2); D2=d2+2; a2=1+A/mu_r*v2
print("flat D2:", sp.simplify(D2-2*(1-th))==0, " flat hull:", sp.simplify(D2*(D2-2*a2)+4*th*(1-th)*(1-A**2)/(1-A**2*r**2))==0)
print("F'' :", sp.simplify(sp.diff(F,r,2)+2*th*(1-A**2)/(1-A*r)**3)==0, " F'(1),F'(-1):", sp.simplify(sp.diff(F,r).subs(r,1)), sp.simplify(sp.diff(F,r).subs(r,-1)))
# Lemma 3.2 numerically: random laminates
rng=np.random.default_rng(0)
err=0
for _ in range(2000):
    Aa=rng.uniform(-.95,.95); rho=rng.uniform(-.99,.99); v=rng.normal(size=2); th_=rng.uniform(0,2*np.pi)
    nu=np.array([np.cos(th_),np.sin(th_)]); mu_=1+Aa*rho; a=np.array([0,1])+Aa/mu_*v
    lp,lm=(1+rho)/2,(1-rho)/2
    delta=-2*np.array([0,1])+2*(a@nu)*nu
    vp,vm=v+lm*delta,v-lp*delta; up,um=vp/(1+Aa),vm/(1-Aa)
    perp=lambda w: np.array([-w[1],w[0]])
    c1=up@nu-um@nu; c2=((1+Aa)*up-(1-Aa)*um)@perp(nu)+2*nu[0]
    u=lp*up+lm*um; m=lp*up-lm*um
    D=2*((1-Aa**2)*m-(rho-Aa)*(u+Aa*m))/(1-rho**2)+2*np.array([0,1])
    err=max(err,abs(c1),abs(c2),np.max(np.abs(D-2*(a@nu)*nu)),np.max(np.abs(u+Aa*m-v)))
print("laminate lemma max error:", err)
