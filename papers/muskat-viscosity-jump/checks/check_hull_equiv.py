# Symbolic proof of the identity behind U_A = O (manuscript, Sec. 3, identification with the relaxed set of [Men22]):
#     D.(D-2a) = 4(1-A^2) g / ((1-th^2)^2 (1+A th)),
# where g = (1-A th)|m-th u|^2 + (1-th^2)(m-th u).(Au+e2) is the function of LM26 (A.1) (Men22, Lemmas 4.1-4.2;
# Szekelyhidi for A=0).  Since the factor is positive for |th|<1, |A|<1, the conditions D.(D-2a)<0 and g<0 agree,
# so U_A is the set O of Men22, written (O) in LM26: a reformulation at fixed Darcy field, not a new set.
import sympy as sp
A,th,u1,u2,m1,m2=sp.symbols('A theta u1 u2 m1 m2',real=True)
u=sp.Matrix([u1,u2]); m=sp.Matrix([m1,m2]); e2=sp.Matrix([0,1])
mu=1+A*th; v=u+A*m
d=2*((1-A**2)*m-(th-A)*v)/(1-th**2); D=d+2*e2; a=e2+(A/mu)*v
mine=sp.together((D.T*(D-2*a))[0])
w=m-th*u; g=(1-A*th)*(w.T*w)[0]+(1-th**2)*(w.T*(A*u+e2))[0]   # Luna-Velasco--Mengual (A.1)
ratio=sp.simplify(mine/g)
factor=4*(1-A**2)/((1-th**2)**2*(1+A*th))
print("D.(D-2a) / g =", sp.factor(ratio))
indep=all(sp.diff(ratio,s)==0 for s in (u1,u2,m1,m2))
print("ratio independent of (u,m):", indep)
assert indep
assert sp.simplify(ratio-factor)==0
print("ratio equals 4(1-A^2)/((1-th^2)^2 (1+A th)):", True)
# positivity for |A|<1, |th|<1: write A=(1-p)/(1+p), th=(1-q)/(1+q) with p,q>0
p,q=sp.symbols('p q',positive=True)
fpq=sp.factor(sp.simplify(factor.subs({A:(1-p)/(1+p),th:(1-q)/(1+q)})))
print("factor with A=(1-p)/(1+p), th=(1-q)/(1+q):", fpq)
assert fpq.is_positive
print("factor is positive for |A|<1, |th|<1:", True)
print("check_hull_equiv.py: all assertions passed")
