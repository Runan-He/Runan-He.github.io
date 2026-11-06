import numpy as np
N=512; h=2*np.pi/N; a=np.arange(N)*h
# test curve z(alpha) = alpha + periodic, complex form
z  = (a + 0.3*np.sin(a)) + 1j*(0.4*np.cos(a)+0.1*np.sin(2*a))
dz = (1 + 0.3*np.cos(a)) + 1j*(-0.4*np.sin(a)+0.2*np.cos(2*a))
d2z= (-0.3*np.sin(a)) + 1j*(-0.4*np.cos(a)-0.4*np.sin(2*a))
L=np.abs(dz); tau=dz/L; n=1j*tau
n2=n.imag; tau2=tau.imag
# kernel of K (adjoint Neumann-Poincare), trapezoidal rule (smooth kernel)
ZA,ZB=np.meshgrid(z,z,indexing='ij'); DA,DB=np.meshgrid(dz,dz,indexing='ij')
with np.errstate(all='ignore'):
    k=(np.abs(DB)/(4*np.pi*np.abs(DA)))*np.imag(DA/np.tan((ZB-ZA)/2))
np.fill_diagonal(k,-np.imag(d2z/dz)/(4*np.pi))
K=k*h
for A in [0.0,0.3,-0.6,0.85]:
    I=np.eye(N)
    sig=np.linalg.solve(I-2*A*K,-2*n2)
    U_sig=-K@sig
    # classical vorticity formulation: (I+2AK)gamma=-2 tau2, U = pv normal comp. of vortex sheet field
    gam=np.linalg.solve(I+2*A*K,-2*tau2)
    U_gam=np.zeros(N)
    for i in range(N):
        j=np.arange((i+1)%2,N,2)  # alternating-point rule for the principal value
        cb=(1/(4j*np.pi))*np.sum(gam[j]*L[j]/np.tan((z[j]-z[i])/2))*2*h
        U_gam[i]=np.real(cb*n[i])
    print(f"A={A:+.2f}  max|U_sig-U_gam|={np.max(np.abs(U_sig-U_gam)):.2e}  int sigma ds={np.sum(sig*L)*h:+.6f} (-4pi={-4*np.pi:+.6f})  "
          f"max|sigma+2(n2+AU)|={np.max(np.abs(sig+2*(n2+A*U_sig))):.1e}  min/max sigma={sig.min():+.3f}/{sig.max():+.3f}")
