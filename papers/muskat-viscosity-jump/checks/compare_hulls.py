import numpy as np
rng=np.random.default_rng(1); e2=np.array([0.,1.]); agree=0; N=200000; ratios=[]
for _ in range(N):
    A=rng.uniform(-.95,.95); th=rng.uniform(-.99,.99); u=rng.normal(size=2)*2; m=rng.normal(size=2)*2
    mu=1+A*th; v=u+A*m
    d=2*((1-A**2)*m-(th-A)*v)/(1-th**2); D=d+2*e2; a=e2+A/mu*v
    mine=D@(D-2*a)            # my set: D.(D-2a)<0
    w=m-th*u; g=(1-A*th)*(w@w)+(1-th**2)*(w@(A*u+e2))   # Luna-Mengual: g<0
    agree+= (mine<0)==(g<0)
    ratios.append(mine/g if abs(g)>1e-9 else np.nan)
    if _<1: pass
r=np.array(ratios); r=r[np.isfinite(r)]
print(f"agreement {agree}/{N}")
# test whether mine = K(theta,A) * g for an explicit factor
A=0.37; th=-0.41
for _ in range(3):
    u=rng.normal(size=2); m=rng.normal(size=2); mu=1+A*th; v=u+A*m
    d=2*((1-A**2)*m-(th-A)*v)/(1-th**2); D=d+2*e2; a=e2+A/mu*v
    w=m-th*u; g=(1-A*th)*(w@w)+(1-th**2)*(w@(A*u+e2))
    print("ratio mine/g =", (D@(D-2*a))/g, " predicted 4(1-A th)/((1-th^2)^2 (1+A th)) =", 4*(1-A*th)/((1-th**2)**2*(1+A*th)))
