"""Regenerate every fitted constant of the final law from named data, in one run.
Data: (1) monodisperse baseline points, (2) polydispersity points, (3) alignment points -- the three lists in fitlaw.py (simulation results divided by 1.015
where marked); (4) thin-limit ripple ratios in ripthin2_G.pkl (paired simulations, flat vs rippled).
Run from the workdir:  python3 calibrate_all.py   -> prints a table and writes calibration_result.json
Fits (all weighted least squares, weights = stated uncertainties):
  B_       : baseline finite-thickness constant, 1 parameter on 12 points
  kappa    : log-normal drift constant, 1 parameter on 8 points (5 thin-limit eigenvalue ratios, 3 finite-aspect simulations)
  (c,gamma): alignment prefactor and domain-size constant, exponent fixed at 2, 2 parameters on 24 points   [revised law]
  gamma0   : domain-size constant with prefactor fixed at 0.625, 1 parameter on 24 points                    [frozen law, for comparison]
  CRIP     : ripple damping constant, form K = 1+(Kell(a)-1)exp(-CRIP a lam), 1 parameter on the thin-limit ripple ratios"""
import sys, io, json, math, pickle, contextlib, numpy as np
from scipy.optimize import least_squares
with contextlib.redirect_stdout(io.StringIO()):
    import fitlaw as f
from mylaw import c0, c1, c2, P
def dil(A, K): return (A+c1+c2/A)/(A/K+c1+c2/A)
out = {}
def fit(name, fun, x0, n):
    r = least_squares(fun, x0); J = r.jac; dof = n-len(x0)
    cov = np.linalg.pinv(J.T@J)*(2*r.cost/dof); out[name] = dict(value=[float(v) for v in r.x], se=[float(v) for v in np.sqrt(np.diag(cov))], chi2=float(2*r.cost), dof=dof)
    print(f"{name:34s} {np.round(r.x,4)} +- {np.round(np.sqrt(np.diag(cov)),4)}  chi2 {2*r.cost:6.1f} / {dof} dof")
FS = 1.015
fit("B_ (baseline)", lambda p: [(f.phiA(A, p[0])/(y/FS)-1)/0.015 for A, y in f.base], [0.5], len(f.base))
fit("kappa (polydispersity)", lambda p: f.res_pol2(p), [0.5], len(f.pol)+len(f.fin))
Kal = lambda S, d, c, g: 1/(1-c*S*(1+g/(d/2)**2))
fit("alignment gamma0, prefactor 0.625", lambda p: [(dil(A, Kal(S, d, 0.625, p[0]))/y-1)/(e/y) for A, d, S, y, e in f.al], [0.3], len(f.al))
fit("alignment (prefactor, gamma) free", lambda p: [(dil(A, Kal(S, d, p[0], p[1]))/y-1)/(e/y) for A, d, S, y, e in f.al], [0.6, 0.3], len(f.al))
r = pickle.load(open('ripthin2_G.pkl', 'rb')); fl = r['flat']; m0 = fl.mean(); e0 = fl.std()/math.sqrt(len(fl)); D = []
for k, v in r.items():
    if k == 'flat': continue
    a, lam = k; R = v.mean()/m0; e = R*math.sqrt((v.std()/math.sqrt(len(v))/v.mean())**2+(e0/m0)**2); D.append((a, lam, R-1, e))
D = np.array(D); ke = np.array([f.Kell(x) for x in D[:, 0]])-1
fit("CRIP (ripple, thin-limit ratios)", lambda p: (ke*np.exp(-p[0]*D[:, 0]*D[:, 1])-D[:, 2])/D[:, 3], [0.6], len(D))
json.dump(out, open("calibration_result.json", "w"), indent=1)
