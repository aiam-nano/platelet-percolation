"""File the law's predictions for Exam B BEFORE any truth simulation is run.
Run from the code/ directory:  python3 exam_B/predict_exam_B.py
Writes exam_B/predictions.json with the predictions, the SHA-256 of spec.json, finallaw3.py, mylaw.py and align.py, and the time of filing,
and exam_B/predictions.sha256 with the SHA-256 of predictions.json.  The hash is to be copied into the project record and, if possible, sent to a
third party (for example an email to a colleague) before the truth is simulated."""
import json, math, hashlib, sys, time, subprocess, numpy as np
sys.path.insert(0, '.')
from align import aligned_normals, shell_P2
import finallaw3 as F
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
spec = json.load(open('exam_B/spec.json'))
KAP = {0.5: 4.733, 0.6: 6.312}; L = 36.0; ell = 6.0
Sloc = {}
for S, kap in KAP.items():
    v = []
    for sd in range(6):
        rng = np.random.default_rng(900 + sd); n = 6000; pos = rng.random((n, 3)) * L
        v.append(shell_P2(pos, aligned_normals(pos, L, ell, kap, rng), L, 2.0)[0])
    Sloc[S] = float(np.mean(v))
pred = {}
for k, c in spec['cases'].items():
    if c['kind'] == 'threshold':
        S = Sloc[c['S_dir']] if c['S_dir'] > 0 else 0.0
        phi = F.phi(A=c['A'], N=c['N'], sig=c['sig'], S=S, d=6.0 if S > 0 else 4.0)
        pred[k] = dict(law_nominal_volume_fraction=phi, law_occupied_fraction=1 - math.exp(-phi), S_used=S)
    else:
        pred[k] = dict(law_ratio=F.Krip(c['a'], c['lam']))
out = dict(filed=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()), hashes=dict(spec=sha('exam_B/spec.json'), finallaw3=sha('finallaw3.py'), mylaw=sha('mylaw.py'), align=sha('align.py')),
           S_local_measured_on_generator=Sloc, predictions=pred)
json.dump(out, open('exam_B/predictions.json', 'w'), indent=1)
open('exam_B/predictions.sha256', 'w').write(sha('exam_B/predictions.json') + '  predictions.json\n')
print(json.dumps(pred, indent=1)); print('predictions.json sha256', sha('exam_B/predictions.json'))
