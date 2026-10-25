"""Score Exam B.  Opens the sealed truth; run only when the freeze review is complete.
usage (from code/):  python3 exam_B/score_exam_B.py
Checks that truth_sealed_B.json matches its recorded hash, that predictions.json matches predictions.sha256 and that the law files match the hashes filed in predictions.json."""
import json, hashlib, math, numpy as np
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
assert sha('exam_B/predictions.json') == open('exam_B/predictions.sha256').read().split()[0], 'predictions changed since filing'
P = json.load(open('exam_B/predictions.json'))
for f, key in (('exam_B/spec.json', 'spec'), ('finallaw3.py', 'finallaw3'), ('mylaw.py', 'mylaw'), ('align.py', 'align')):
    ok = sha(f) == P['hashes'][key]; print(f'{f:22s} unchanged since filing: {ok}')
assert sha('exam_B/truth_sealed_B.json') == open('exam_B/truth_sealed_B.sha256').read().split()[0], 'truth file changed'
T = json.load(open('exam_B/truth_sealed_B.json')); errs = []
for k, p in P['predictions'].items():
    if 'law_occupied_fraction' in p: pv = p['law_occupied_fraction']
    else: pv = p['law_ratio']
    e = 100 * (pv / T[k] - 1); errs.append(abs(e)); print(f'{k:4s} law {pv:.5f} truth {T[k]:.5f} error {e:+.1f}%')
print(f'mean absolute error over {len(errs)} cases: {np.mean(errs):.2f}%')
