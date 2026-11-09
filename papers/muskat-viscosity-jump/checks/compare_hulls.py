# Numerical comparison of the relaxed set U_A = {|rho|<1, D.(D-2a)<0} of the manuscript (Sec. 3) with the set O of
# Mengual [Men22], written (O) in LM26 = Luna-Velasco & Mengual, arXiv:2610.34504, through the function
# g = (1-A th)|m-th u|^2 + (1-th^2)(m-th u).(Au+e2) of LM26 (A.1).  The two sets coincide because
#     D.(D-2a) = 4(1-A^2) g / ((1-th^2)^2 (1+A th)),
# an exact identity (proved symbolically in check_hull_equiv.py); this is a reformulation of O, not a new set.
import numpy as np
rng=np.random.default_rng(1); e2=np.array([0.,1.]); agree=0; N=200000

def states(A,th,u,m):
    mu=1+A*th; v=u+A*m
    d=2*((1-A**2)*m-(th-A)*v)/(1-th**2); D=d+2*e2; a=e2+A/mu*v
    mine=D@(D-2*a)            # manuscript: D.(D-2a)<0
    w=m-th*u; g=(1-A*th)*(w@w)+(1-th**2)*(w@(A*u+e2))   # Men22 / LM26 (A.1): g<0
    return mine,g

def factor(A,th):             # predicted ratio D.(D-2a)/g
    return 4*(1-A**2)/((1-th**2)**2*(1+A*th))

maxrel=0.0
for _ in range(N):
    A=rng.uniform(-.95,.95); th=rng.uniform(-.99,.99); u=rng.normal(size=2)*2; m=rng.normal(size=2)*2
    mine,g=states(A,th,u,m)
    agree+= (mine<0)==(g<0)
    if abs(g)>1e-6:
        maxrel=max(maxrel,abs(mine/g-factor(A,th))/factor(A,th))
print(f"agreement of the two sign conditions: {agree}/{N}")
print(f"max relative deviation of D.(D-2a)/g from 4(1-A^2)/((1-th^2)^2 (1+A th)) over the sample: {maxrel:.1e}")
assert agree==N
assert maxrel<1e-8

# the ratio at a fixed (A, th), for a few random (u, m)
A=0.37; th=-0.41
for _ in range(3):
    u=rng.normal(size=2); m=rng.normal(size=2)
    mine,g=states(A,th,u,m)
    print("ratio D.(D-2a)/g =", mine/g, "  predicted 4(1-A^2)/((1-th^2)^2 (1+A th)) =", factor(A,th))
    assert abs(mine/g-factor(A,th))<1e-9*factor(A,th)
assert factor(A,th)>0
print("compare_hulls.py: all assertions passed")
