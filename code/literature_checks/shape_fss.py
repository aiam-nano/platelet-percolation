"""First-wrap threshold of zero-thickness plates (disc, square, equilateral triangle) with the same estimator,
box sizes L (in units of the equal-area radius r, so eta = (4/3) pi N r^3 / L^3 as in Yi and Tawerghi).
usage: python3 shape_fss.py shape L nreal seed [Nfac]"""
import sys, time, math, json
import numpy as np
from sim2e import *
from poly_test import polygon

def run(shape, L, nreal, seed, Nfac=0.3):
    rng = np.random.default_rng(seed); out = []
    if shape == 'disc':
        mode, V, br = 0, np.zeros((1, 3)), 1.0
    else:
        mode, V = 4, polygon(shape); br = float(np.linalg.norm(V, axis=1).max())
    for _ in range(nreal):
        n = int(Nfac * L ** 3)
        pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng)
        rad = np.ones(n); hth = np.zeros(n); brad = np.full(n, br)
        Nc = run_once(mode, pos, rot, rad, hth, L, brad) if mode == 0 else \
            first_perc(mode, pos, rot, rad, hth, L, max(3, int(L // (2 * br))), brad, V, EMPTY['E'], 0, 0., 0., 0., 1., EMPTY['xtab'], EMPTY['utab'])
        if Nc > 0: out.append(Nc * 4 / 3 * math.pi / L ** 3)
    return np.array(out)

if __name__ == "__main__":
    shape, L, nreal, seed = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    t0 = time.time(); e = run(shape, int(L), nreal, seed)
    row = dict(shape=shape, L=int(L), n=len(e), eta=float(e.mean()), err=float(e.std(ddof=1) / math.sqrt(len(e))), seed=seed, secs=round(time.time() - t0))
    print(json.dumps(row), flush=True)
    open('shape_fss.jsonl', 'a').write(json.dumps(row) + '\n')
