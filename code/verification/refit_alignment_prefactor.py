"""Provenance check: is the alignment prefactor 0.625 supported by the alignment data it is used on, or only assumed?
Refits the alignment ratios (24 simulated points as listed in fitlaw.py) with the prefactor free, exponent fixed at 2 (as in finallaw2)."""
import sys, math, numpy as np
sys.path.insert(0, "..")
from scipy.optimize import least_squares
from mylaw import c0, c1, c2
al=[(33.333,2,.119,1.0995,.012),(33.333,2,.146,1.1358,.012),(33.333,2,.166,1.141,.012),(33.333,2,.193,1.177,.012),
    (33.333,4,.282,1.206,.012),(33.333,4,.314,1.2124,.012),
    (33.333,6,.180,1.0858,.012),(33.333,6,.268,1.1761,.012),(33.333,6,.355,1.236,.012),
    (33.333,8,.252,1.1431,.012),(33.333,8,.398,1.268,.012),(33.333,8,.478,1.3736,.012),
    (30,4,.239,1.1672,.012),(30,4,.283,1.1932,.012),(30,8,.339,1.2216,.012),(30,8,.397,1.263,.012),(30,8,.476,1.3587,.012),
    (10,4,.181,1.083,.012),(10,4,.236,1.1076,.012),
    (1e6,3,.044,1.0425,.01),(1e6,3,.217,1.188,.01),(1e6,3,.330,1.361,.01),(1e6,6,.073,1.0588,.01),(1e6,6,.356,1.276,.01)]
def dil(A,K): return (A+c1+c2/A)/(A/K+c1+c2/A)
def Kal(S,d,c,g,pw=2.0): return 1/(1-c*S*(1+g/(d/2)**pw))
def chi2(c,g,pw=2.0):
    r=[(dil(A,Kal(S,d,c,g,pw))/y-1)/(e/y) for A,d,S,y,e in al]; return np.array(r)
for label,fun,x0,fixed in [
    ("c fixed 0.625, g free",lambda p: chi2(0.625,p[0]),[0.3],None),
    ("c free, g free",lambda p: chi2(p[0],p[1]),[0.6,0.3],None),
    ("c free, g free, exponent free",lambda p: chi2(p[0],p[1],p[2]),[0.6,0.3,2.0],None)]:
    r=least_squares(fun,x0); J=r.jac; cov=np.linalg.pinv(J.T@J)*(2*r.cost/(len(al)-len(x0)))
    print(f"{label}: params {np.round(r.x,3)} +- {np.round(np.sqrt(np.diag(cov)),3)}  chi2={2*r.cost:.1f} on {len(al)-len(x0)} dof")
