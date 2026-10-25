"""Extends fss_thin.py to larger boxes (nfac reduced because the threshold needs only about 0.23 L^3 discs).
Usage: python3 fss_thin_large.py L nreal seed out.json"""
import sys, json, math, time, numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "..")
from sim2 import random_rotations, run_once
L = float(sys.argv[1]); nreal = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3])); v = []; t0 = time.time()
for _ in range(nreal):
    n = int(0.30 * L**3); pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng)
    Nc = run_once(0, pos, rot, np.ones(n), np.zeros(n), L, np.ones(n))
    if Nc > 0: v.append(Nc / L**3 * 4 * math.pi / 3)
v = np.array(v); json.dump(v.tolist(), open(sys.argv[4], "w"))
print(f"L={L:4.0f} n={len(v)} eta={v.mean():.4f} +- {v.std()/math.sqrt(len(v)):.4f} sd={v.std():.4f} [{time.time()-t0:.0f}s]", flush=True)
