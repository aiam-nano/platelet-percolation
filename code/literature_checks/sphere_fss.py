"""Finite-size extrapolation of the sphere threshold with the winding-number estimator (control for sim2e.py mode 3, aspect ratio 1)."""
import numpy as np, math
from ell_test import run
res=[]
for L,nr in ((10,200),(14,120),(20,80),(28,40),(40,20),(56,10)):
    nc=run(3,1.0,L,nr,777+L,0.5); eta=nc*4/3*math.pi
    res.append((L,eta.mean(),eta.std()/math.sqrt(len(eta)))); print(L,round(eta.mean(),4),round(res[-1][2],4),flush=True)
nu=0.8765
x=np.array([r[0]**(-1/nu) for r in res]); y=np.array([r[1] for r in res]); w=1/np.array([r[2] for r in res])**2
A=np.vstack([np.ones_like(x),x]).T; beta=np.linalg.solve(A.T@np.diag(w)@A,A.T@np.diag(w)@y)
print("extrapolated nominal eta_c =",round(beta[0],4),"(accepted 0.3418)")
