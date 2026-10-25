"""Per-axis winding and estimator check for thin discs: N_any versus the first wrap along each axis. Run from workdir: python3 verification/axes_check.py L nreal seed out.json"""
import sys, json, math, numpy as np
sys.path.insert(0, ".")
from sim2 import random_rotations
from sim2_axes import first_perc_axes
L = float(sys.argv[1]); nreal = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3])); rows = []
for _ in range(nreal):
    n = int(0.40*L**3); pos = rng.random((n, 3))*L; rot = random_rotations(n, rng)
    r = first_perc_axes(0, pos, rot, np.ones(n), np.zeros(n), L, max(3, int(L//2)), np.ones(n)); rows.append(r.tolist())
a = np.array(rows, float); eta = a/L**3*4*math.pi/3
json.dump(rows, open(sys.argv[4], "w"))
names = ["any", "x", "y", "z"]
for k in range(4):
    print(f"L={L:g} N_{names[k]:3s}: eta mean {eta[:,k].mean():.4f} +- {eta[:,k].std()/math.sqrt(len(a)):.4f}  median {np.median(eta[:,k]):.4f}")
ax = eta[:, 1:].mean(1); print(f"mean over axes of N_axis: {ax.mean():.4f}; max over axes {eta[:,1:].max(1).mean():.4f}")
