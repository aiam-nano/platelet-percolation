"""Leave-one-point-out cross-validation of each single-constant fit, using exactly the weighted residuals of calibrate_all.py.
Reports RMS of the normalised residual (data-weighted) and RMS percentage error in-sample and for the held-out point. Run from the workdir."""
import sys, io, math, pickle, contextlib, numpy as np
sys.path.insert(0, ".")
from scipy.optimize import least_squares
with contextlib.redirect_stdout(io.StringIO()):
    import fitlaw as f
from mylaw import c0, c1, c2, P
FS = 1.015
def dil(A, K): return (A+c1+c2/A)/(A/K+c1+c2/A)
Kal = lambda S, d, c, g: 1/(1-c*S*(1+g/(d/2)**2))
r = pickle.load(open('ripthin2_G.pkl', 'rb')); fl = r['flat']; m0 = fl.mean(); e0 = fl.std()/math.sqrt(len(fl)); D = []
for k, v in r.items():
    if k == 'flat': continue
    a, lam = k; R = v.mean()/m0; e = R*math.sqrt((v.std()/math.sqrt(len(v))/v.mean())**2+(e0/m0)**2); D.append((a, lam, R-1, e))
fam = {  # name: (data, x0, (weighted residual, percentage error))
 "base B_": (f.base, [0.5], lambda p, x: (f.phiA(x[0], p[0])/(x[1]/FS)-1)/0.015, lambda p, x: 100*(f.phiA(x[0], p[0])/(x[1]/FS)-1)),
 "alignment gamma (prefactor 0.625)": (f.al, [0.3], lambda p, x: (dil(x[0], Kal(x[2], x[1], 0.625, p[0]))/x[3]-1)/(x[4]/x[3]), lambda p, x: 100*(dil(x[0], Kal(x[2], x[1], 0.625, p[0]))/x[3]-1)),
 "alignment (prefactor, gamma) free": (f.al, [0.6, 0.3], lambda p, x: (dil(x[0], Kal(x[2], x[1], p[0], p[1]))/x[3]-1)/(x[4]/x[3]), lambda p, x: 100*(dil(x[0], Kal(x[2], x[1], p[0], p[1]))/x[3]-1)),
 "ripple CRIP (ratio-1, absolute)": (D, [0.6], lambda p, x: (( f.Kell(x[0])-1)*math.exp(-p[0]*x[0]*x[1])-x[2])/x[3], lambda p, x: 100*(((f.Kell(x[0])-1)*math.exp(-p[0]*x[0]*x[1])-x[2])/(1+x[2]))),
 "polydispersity kappa": ([("thin",)+tuple(x) for x in f.pol]+[("fin",)+tuple(x) for x in f.fin], [0.5],
     lambda p, x: ((1+p[0]*x[2]**2)/x[3]-1)/(x[4]/x[3]) if x[0]=="thin" else (P(x[2], x[1], p[0])/x[3]-1)/0.02,
     lambda p, x: 100*((1+p[0]*x[2]**2)/x[3]-1) if x[0]=="thin" else 100*(P(x[2], x[1], p[0])/x[3]-1)),
}
for name, (data, x0, wres, pct) in fam.items():
    data = list(data); n = len(data)
    full = least_squares(lambda p: [wres(p, x) for x in data], x0)
    ins = math.sqrt(np.mean([wres(full.x, x)**2 for x in data])); insp = math.sqrt(np.mean([pct(full.x, x)**2 for x in data]))
    lo, lop = [], []
    for i in range(n):
        rr = least_squares(lambda p: [wres(p, x) for j, x in enumerate(data) if j != i], x0); lo.append(wres(rr.x, data[i])); lop.append(pct(rr.x, data[i]))
    print(f"{name:36s} n={n:2d} value={np.round(full.x,3)}  in-sample RMS {ins:.2f} sigma / {insp:.2f}%   leave-one-out RMS {math.sqrt(np.mean(np.square(lo))):.2f} sigma / {math.sqrt(np.mean(np.square(lop))):.2f}%  worst {max(abs(np.array(lop))):.1f}%")
