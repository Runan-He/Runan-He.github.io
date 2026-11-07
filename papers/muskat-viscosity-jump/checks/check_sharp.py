import sympy as sp
A,rm,k,s,g=sp.symbols('A rho_m kappa sigma0 gamma0',real=True)
mu=1+A*rm
Dn=-s-2*(1-A*rm)*k
an=-s/2-A**2*(1-rm**2)*k/mu
at=-(1-A**2)*g/(2*mu)
expr=Dn*(Dn-2*an)-at**2            # D.(D-2a) with D_tau=a_tau
target=(2*(1-A**2)/mu)*(k*(2*(1-A*rm)*k+s)-(1-A**2)*g**2/(8*mu))
print("sharp-range identity:", sp.simplify(expr-target)==0)
# optimality of the tilt D_tau=a_tau: min over D_tau of D_tau(D_tau-2a_tau) is -a_tau^2
Dt=sp.symbols('D_tau',real=True)
print("tilt completes the square:", sp.simplify(Dt*(Dt-2*at)+at**2-(Dt-at)**2)==0)
kmax=(-s+sp.sqrt(s**2+(1-A*rm)*(1-A**2)*g**2/mu))/(4*(1-A*rm))
print("kappa_max is the root:", sp.simplify((kmax*(2*(1-A*rm)*kmax+s)-(1-A**2)*g**2/(8*mu)))==0)
print("gamma0=0 reduces to old range:", sp.simplify(kmax.subs(g,0).subs(s,-sp.Symbol('S',positive=True)))) 
