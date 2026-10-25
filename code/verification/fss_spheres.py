"""Control experiment: the same percolation pipeline (cell list, union-find with periodic winding) applied to overlapping SPHERES, whose
infinite-system threshold is well established: phi_c = 0.2895 +- 0.0005 (Lorenz and Ziff 2001; Rintoul and Torquato 1997), i.e. eta_c = n_c (4 pi/3) R^3 = -ln(1-0.2895) = 0.3418.
Usage: python3 fss_spheres.py nreal seed out.json"""
import sys, json, math, time, numpy as np
from sim2_sphere_check import random_rotations, run_once
nreal = int(sys.argv[1]); rng = np.random.default_rng(int(sys.argv[2])); res = {}
for L in (6.0, 8.0, 12.0, 16.0, 24.0, 32.0):
    v = []; t0 = time.time()
    for _ in range(nreal):
        n = int(0.9 * L**3); pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng)
        Nc = run_once(3, pos, rot, np.ones(n), np.zeros(n), L, np.ones(n))
        if Nc > 0: v.append(Nc / L**3 * 4 * math.pi / 3)
    v = np.array(v); res[str(L)] = v.tolist()
    print(f"L={L:4.0f} n={len(v)} eta={v.mean():.4f} +- {v.std()/math.sqrt(len(v)):.4f} [{time.time()-t0:.0f}s]", flush=True)
    json.dump(res, open(sys.argv[3], "w"))
