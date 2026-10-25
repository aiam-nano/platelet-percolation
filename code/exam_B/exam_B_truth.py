"""Truth simulation for Exam B.  Run only AFTER predictions.json has been filed and its hash recorded.
usage (from code/):  python3 exam_B/exam_B_truth.py  [case ...]
Writes exam_B/truth_B_raw.pkl and exam_B/truth_sealed_B.json (values rounded to five figures) with its SHA-256 in truth_sealed_B.sha256.
It prints only the case names and timings, never the values, so that the operator opens the file at scoring time."""
import sys, json, math, time, pickle, hashlib, numpy as np
sys.path.insert(0, '.')
from sim2 import run_once, random_rotations, rotation_with_normal
from align import aligned_normals
import gap_rip
spec = json.load(open('exam_B/spec.json'))
KAP = {0.5: 4.733, 0.6: 6.312}
def lognormal_R(n, s, rng, rcap=6.0):
    R = rng.lognormal(-0.5 * s * s, s, size=n); return R[R < rcap + 1e-12]
def one(c, L, nfac, rng):
    n = int(nfac * L ** 3)
    R = lognormal_R(n, c['sig'], rng) if c['sig'] > 0 else np.ones(n); n = len(R)
    t = np.full(n, c['N'] * 2.0 / c['A'])
    pos = rng.random((n, 3)) * L
    if c['S_dir'] > 0:
        nrm = aligned_normals(pos, L, 6.0, KAP[c['S_dir']], rng); rot = rotation_with_normal(nrm, rng)
    else:
        rot = random_rotations(n, rng)
    hth = t / 2; brad = np.sqrt(R ** 2 + hth ** 2)
    Nc = run_once(1, pos, rot, R, hth, L, brad)
    if Nc < 0: return None
    return Nc / L ** 3 * np.mean(math.pi * R[:Nc] ** 2 * t[:Nc])
if __name__ == '__main__':
    cases = sys.argv[1:] or list(spec['cases']); res = {}; t0 = time.time()
    for k in cases:
        c = spec['cases'][k]; seed = 71000 + sum(map(ord, k)); rng = np.random.default_rng(seed)
        if c['kind'] == 'threshold':
            v = [x for x in (one(c, 32.0, 2.5, rng) for _ in range(48)) if x is not None]
            res[k] = np.array(v)
        else:
            gap_rip.RCAP = 6.0; o = gap_rip.paired(32.0, c['sig'], c['a'], c['lam'], 48, seed); res[k] = o
        pickle.dump(res, open('exam_B/truth_B_raw.pkl', 'wb')); print(k, 'done', f'[{time.time()-t0:.0f}s]', flush=True)
    sealed = {}
    for k, a in res.items():
        c = spec['cases'][k]
        if c['kind'] == 'threshold':
            x = a.mean() / 1.015; sealed[k] = round(1 - math.exp(-x), 5)
        else:
            sealed[k] = round(float(a[:, 1].mean() / a[:, 0].mean()), 4)
    raw = json.dumps(sealed, sort_keys=True).encode()
    open('exam_B/truth_sealed_B.json', 'wb').write(raw)
    open('exam_B/truth_sealed_B.sha256', 'w').write(hashlib.sha256(raw).hexdigest() + '  truth_sealed_B.json\n')
    print('sealed; sha256', hashlib.sha256(raw).hexdigest())
