import numpy as np
N=512; h=2*np.pi/N; a=np.arange(N)*h
def curve(name):
    if name=='turned':   # overturned periodic interface: z = (alpha - 1.5 sin, 0.8 sin)
        z=(a-1.5*np.sin(a))+1j*(0.8*np.sin(a)); dz=(1-1.5*np.cos(a))+1j*(0.8*np.cos(a)); d2z=(1.5*np.sin(a))+1j*(-0.8*np.sin(a))
    else:                # graph from check_sigma.py
        z=(a+0.3*np.sin(a))+1j*(0.4*np.cos(a)+0.1*np.sin(2*a)); dz=(1+0.3*np.cos(a))+1j*(-0.4*np.sin(a)+0.2*np.cos(2*a)); d2z=(-0.3*np.sin(a))+1j*(-0.4*np.cos(a)-0.4*np.sin(2*a))
    return z,dz,d2z
def ops(z,dz,d2z):
    L=np.abs(dz); tau=dz/L; n=1j*tau
    ZA,ZB=np.meshgrid(z,z,indexing='ij'); DA,DB=np.meshgrid(dz,dz,indexing='ij')
    with np.errstate(all='ignore'):
        k=(np.abs(DB)/(4*np.pi*np.abs(DA)))*np.imag(DA/np.tan((ZB-ZA)/2))
    np.fill_diagonal(k,-np.imag(d2z/dz)/(4*np.pi))
    K=k*h
    def T(phi):   # tangential pv component of the source-sheet field  <C[z;phi n]>.tau
        out=np.zeros(N)
        for i in range(N):
            j=np.arange((i+1)%2,N,2)
            cb=-(1/(4*np.pi))*np.sum(phi[j]*L[j]/np.tan((z[j]-z[i])/2))*2*h
            out[i]=np.real(cb*tau[i])
        return out
    return L,tau,n,K,T
for name in ['graph','turned']:
    z,dz,d2z=curve(name); L,tau,n,K,T=ops(z,dz,d2z); I=np.eye(N)
    for A in [0.0,0.3,-0.6,0.9]:
        sig=np.linalg.solve(I-2*A*K,-2*n.imag)
        gam=np.linalg.solve(I+2*A*K,-2*tau.imag)
        lhs=tau.imag-A*T(sig); rhs=-(1-A**2)*gam/2
        st=sig>=0
        pinch_margin=np.min(np.abs(gam[st])) if st.any() else np.inf
        rm=A; mu=1+A*rm
        kmax=(-sig+np.sqrt(sig**2+(1-A*rm)*(1-A**2)*gam**2/mu))/(4*(1-A*rm))
        print(f"{name:6s} A={A:+.1f} | max|tau2-A*T(sig)+(1-A^2)gam/2|={np.max(np.abs(lhs-rhs)):.1e} | "
              f"stable fraction={st.mean():.2f} | min|gamma| on stable set={pinch_margin:.3f} | "
              f"min kappa_max={kmax.min():.3f} (old range min={np.min(np.maximum(-sig,0))/(2*(1-A*rm)):.3f})")
