"""Finite-size scaling of the first-wrapping-cluster threshold for penetrable zero-thickness discs.
Reference infinite-system value eta_c = 0.9614 +- 0.0005 (Yi and Tawerghi 2009).  eta = (4 pi/3) R^3 n_c.
Usage: python3 fss_thin.py nreal seed out.json"""
import sys, json, math, time, numpy as np
sys.path.insert(0, "..")
from sim2 import random_rotations, run_once
nreal = int(sys.argv[1]); seed = int(sys.argv[2]); out = sys.argv[3]
rng = np.random.default_rng(seed); res = {}
for L in (8.0, 12.0, 16.0, 24.0, 32.0, 48.0):
    t0 = time.time(); v = []
    for _ in range(nreal):
        nfac = 1.5
        n = int(nfac * L**3); pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng)
        Nc = run_once(0, pos, rot, np.ones(n), np.zeros(n), L, np.ones(n))
        if Nc > 0: v.append(Nc / L**3 * 4 * math.pi / 3)
    v = np.array(v); res[str(L)] = v.tolist()
    print(f"L={L:4.0f} n={len(v)} eta={v.mean():.4f} +- {v.std()/math.sqrt(len(v)):.4f} sd={v.std():.4f} [{time.time()-t0:.0f}s]", flush=True)
    json.dump(res, open(out, "w"))
