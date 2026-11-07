import sympy as sp
A,th,u1,u2,m1,m2=sp.symbols('A theta u1 u2 m1 m2',real=True)
u=sp.Matrix([u1,u2]); m=sp.Matrix([m1,m2]); e2=sp.Matrix([0,1])
mu=1+A*th; v=u+A*m
d=2*((1-A**2)*m-(th-A)*v)/(1-th**2); D=d+2*e2; a=e2+(A/mu)*v
mine=sp.together((D.T*(D-2*a))[0])
w=m-th*u; g=(1-A*th)*(w.T*w)[0]+(1-th**2)*(w.T*(A*u+e2))[0]   # Luna-Velasco--Mengual
ratio=sp.simplify(mine/g)
print("D.(D-2a) / g =", sp.factor(ratio))
print("ratio independent of (u,m):", all(sp.diff(ratio,s)==0 for s in (u1,u2,m1,m2)))
