"""Prolate spheroids, aspect ratio eps (polar/equatorial), larger boxes than prolate_test.py. usage: python3 prolate_fss.py eps L nreal seed"""
import sys, math, json, time
import numpy as np
from sim2e import *
eps, L, nreal, seed = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed); out = []; t0 = time.time(); vol = 4 / 3 * math.pi * eps
for _ in range(nreal):
    n = int(0.40 / vol * L ** 3)
    pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng)
    rad = np.ones(n); hth = np.full(n, eps); brad = np.full(n, eps)
    Nc = run_once(3, pos, rot, rad, hth, L, brad)
    if Nc > 0: out.append(Nc / L ** 3 * vol)
ph = np.array(out); occ = 1 - np.exp(-ph)
row = dict(eps=eps, L=L, n=len(ph), phi_nominal=float(ph.mean()), occupied=float(occ.mean()), err=float(occ.std(ddof=1) / math.sqrt(len(ph))), seed=seed, secs=round(time.time() - t0))
print(json.dumps(row), flush=True); open('prolate_fss.jsonl', 'a').write(json.dumps(row) + '\n')
