"""T10: cell-list correctness in boxes where the grid has more than three cells per side. With ncell = 3 the 27-cell neighbourhood covers the whole periodic box,
so it acts as an all-pairs oracle; the default (finer) grid must give the identical first-percolation index on identical configurations.
Run from the workdir (or from verification/): python3 verification/test_larger_box.py"""
import sys, os, math, numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from sim2 import random_rotations, run_once
def lognormal(n, s, rng, cap=6.0):
    R = rng.lognormal(-0.5*s*s, s, size=n); return R[R < cap]
bad = 0; total = 0
for name, L, sig, mode, tfull in [("thin discs, monodisperse", 24.0, 0.0, 0, 0.0), ("cylinders A=20, monodisperse", 20.0, 0.0, 1, 0.1), ("log-normal 0.3, thin cylinders", 32.0, 0.3, 1, 0.04), ("log-normal 0.5, thin cylinders", 40.0, 0.5, 1, 0.04)]:
    rng = np.random.default_rng(4242); res = []
    for k in range(4):
        n = int((0.45 if mode == 0 else 1.2)*L**3/ (1 if sig == 0 else 1)); R = np.ones(n) if sig == 0 else lognormal(n, sig, rng); n = len(R)
        pos = rng.random((n, 3))*L; rot = random_rotations(n, rng); hth = np.full(n, tfull/2); brad = np.sqrt(R**2+hth**2)
        a = run_once(mode, pos, rot, R, hth, L, brad, ncell=3); b = run_once(mode, pos, rot, R, hth, L, brad)
        res.append((a, b)); total += 1; bad += (a != b)
    nc = max(3, int(L // (2*brad.max())))
    print(f"{name:34s} L={L:g} default ncell={nc:3d}: ncell=3 vs default first-percolation index {res}  {'OK' if all(a==b for a,b in res) else 'MISMATCH'}")
print("T10", "PASS" if bad == 0 else f"FAIL ({bad} of {total} differ)")
