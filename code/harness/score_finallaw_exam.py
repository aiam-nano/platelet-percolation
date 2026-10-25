"""Score finallaw2 on the sealed 12-case exam. Reproduces 06_FinalLaw_Exam_Scoring.md.
Run from the workdir that holds finallaw2.py, mylaw.py, align.py, sim2.py and exam/truth_sealed.json.
Disclosed post hoc choices: S_local measured on the generator (shell r<2R), d = 6 (generator domain size), nominal -> occupied conversion."""
import json, math, hashlib, numpy as np
from align import aligned_normals, shell_P2
import finallaw2 as F
raw = open('exam/truth_sealed.json', 'rb').read()
print('truth sha256', hashlib.sha256(raw).hexdigest())
t = json.loads(raw)
KAP = {0.5: 4.733, 0.6: 6.312}; L = 36.0; ell = 6.0
Sloc = {}
for S, kap in KAP.items():
    v = []
    for sd in range(6):
        rng = np.random.default_rng(900+sd); n = 6000; pos = rng.random((n, 3))*L
        v.append(shell_P2(pos, aligned_normals(pos, L, ell, kap, rng), L, 2.0)[0])
    Sloc[S] = float(np.mean(v)); print(f'S_dir {S}: S_loc {Sloc[S]:.4f} (ratio to S_dir^2 {Sloc[S]/S**2:.3f})')
cases = {'F1': dict(A=50), 'F2': dict(A=8), 'F3': dict(A=60, N=4), 'F4': dict(A=80, sig=.4),
         'F8': dict(A=40, S=Sloc[0.6], d=6), 'F9': dict(A=60, N=3, S=Sloc[0.6], d=6),
         'F12': dict(A=60, N=3, sig=.3, S=Sloc[0.5], d=6)}
errs = []
for k, a in cases.items():
    occ = 1-math.exp(-F.phi(**a)); e = 100*(occ/t[k]-1); errs.append(abs(e))
    print(f'{k:4s} law {occ:.5f} truth {t[k]:.5f} err {e:+.1f}%')
for k in ('F10', 'F11'):   # rippled / flat threshold ratio at fixed sheet (material) area, a=1.2, lambda=1.0
    r = F.Krip(1.2, 1.0); e = 100*(r/t[k]-1); errs.append(abs(e)); print(f'{k:4s} law ratio {r:.4f} truth {t[k]:.4f} err {e:+.1f}%')
print(f'MAE over {len(errs)} covered cases: {np.mean(errs):.2f}%   (F5, F6, F7 outside the law)')
