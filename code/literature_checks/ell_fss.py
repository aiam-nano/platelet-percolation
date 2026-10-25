"""Finite-size extrapolation of the oblate-spheroid thin limit with the first-wrap estimator (mode 3 of sim2e.py).
Boxes in units of the equatorial radius; eta = (4/3) pi n R^3.   usage: python3 ell_fss.py u L nreal seed"""
import sys, math, json, time
import numpy as np
from ell_test import run
u, L, nreal, seed = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
t0 = time.time(); nc = run(3, u, L, nreal, seed, 0.3); eta = nc * 4 / 3 * math.pi
row = dict(u=u, L=L, n=len(eta), eta=float(eta.mean()), err=float(eta.std(ddof=1) / math.sqrt(len(eta))), seed=seed, secs=round(time.time() - t0))
print(json.dumps(row), flush=True); open('ell_fss.jsonl', 'a').write(json.dumps(row) + '\n')
