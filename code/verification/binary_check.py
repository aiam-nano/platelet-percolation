"""Published-value check: Yi & Tawerghi (2009) report that a 50:50 binary mixture of thin discs (radius ratio lambda) percolates at
essentially the same eta (defined with <R^3>) as equal discs (<1% difference for lambda ~ 0.5). Test with sim2 at L=40."""
import sys, json, math, numpy as np
sys.path.insert(0, "..")
from sim2 import random_rotations, run_once
L = 40.0; nreal = 80; out = {}
for lam in (1.0, 0.5, 0.25):
    rng = np.random.default_rng(311)
    r0 = 1.0
    r = (1 + lam) * r0 / (4 * (1 + lam**3)) ** (1/3)
    R1, R2 = 2*lam*r/(1+lam), 2*r/(1+lam)
    vals = []
    for _ in range(nreal):
        n = int(0.40 * L**3)
        pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng)
        rad = np.where(np.arange(n) % 2 == 0, R1, R2).astype(float)
        rng.shuffle(rad)
        Nc = run_once(0, pos, rot, rad, np.zeros(n), L, rad)
        if Nc > 0:
            vals.append(Nc / L**3 * 4*math.pi/3 * np.mean(rad[:Nc]**3))
    v = np.array(vals); out[str(lam)] = (v.mean(), v.std()/math.sqrt(len(v)), len(v))
    print(f"lambda={lam}: eta(<R^3>)={v.mean():.4f} +- {v.std()/math.sqrt(len(v)):.4f} (n={len(v)})", flush=True)
json.dump(out, open("binary_check.json", "w"))
base = out["1.0"][0]
for k, (m, e, n) in out.items():
    print(f"  lambda={k}: relative difference to equal discs = {100*(m/base-1):+.2f}%")
