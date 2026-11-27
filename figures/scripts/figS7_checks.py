from fig_common import *
import time
from sim4 import run_once
from sim2 import random_rotations
# (a) first-wrap count distributions
def runs(kind, L, nrep, A=30.0):
    out = []
    for s in range(nrep):
        rng = np.random.default_rng(5000 + s); n = int(2.5 * L ** 3); pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng); rad = np.ones(n)
        if kind == 'disc': Nc = run_once(0, pos, rot, rad, np.zeros(n), L, rad)
        else: hth = np.full(n, 1.0 / A * 1.0); brad = np.sqrt(rad ** 2 + hth ** 2); Nc = run_once(1, pos, rot, rad, hth, L, brad)
        if Nc > 0: out.append(Nc)
    return np.array(out)
t0 = time.time(); R = {}
for L in (8, 16, 32):
    R[L] = runs('disc', L, 200 if L < 32 else 60); print(L, len(R[L]), R[L].mean() / L ** 3, R[L].std() / R[L].mean(), time.time() - t0, flush=True)
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw=dict(wspace=0.55))
ax = axs[0]; cl = {8: '#bdbdbd', 16: '#4d4d4d', 32: '#000000'}
for L in (8, 16, 32):
    nd = R[L] / L ** 3; ax.hist(nd, bins=np.linspace(0.12, 0.36, 30), density=True, histtype='step', color=cl[L], lw=1.1, label=f'$L={L}$')
ax.set_xlabel('first-wrap density $n_c R^3$'); ax.set_ylabel('density'); ax.set_xlim(0.1, 0.5); ax.legend(frameon=False, loc='upper right', handlelength=1.2)
write_csv('figS7a_wrap_counts', ['L', 'index', 'first_wrap_count'], [(L, i, int(v)) for L in R for i, v in enumerate(R[L])])
# (b) closed-form parity
ref = [('discs', 9.8696, 9.8715, 0.0024), ('720-gons', 9.8696, 9.870, 0.006), ('squares', 11.137, 11.123, 0.008), ('triangles', 12.693, 12.696, 0.012),
       ('sphere\n$c/a=1$', 33.5103, 33.5148, 0.0146), ('spheroid\n$c/a=1/3$', 15.2297, 15.2320, 0.0099), ('spheroid\n$c/a=0.1$', 11.0511, 11.0408, 0.0084), ('spheroid\n$c/a=0.01$', 9.9591, 9.9543, 0.0080)]
ax = axs[1]; y = np.arange(len(ref))[::-1]
ax.errorbar([r[2] / r[1] - 1 for r in ref], y, xerr=[r[3] / r[1] for r in ref], fmt='o', ms=4.5, color=COL['sim'], mfc='white', mew=1.0, lw=1.0)
ax.axvline(0, color=COL['base'], lw=1.0); ax.set_yticks(y); ax.set_yticklabels([r[0] for r in ref]); ax.set_xlabel('Fractional difference from closed form'); ax.set_xlim(-0.004, 0.004)
ax.set_xticks([-0.003, 0, 0.003]); ax.set_xticklabels(['$-0.3$%', '0', '$+0.3$%'])
write_csv('figS7b_closed_forms', ['case', 'closed_form', 'simulated', 'se'], [(r[0].replace('\n', ' '), r[1], r[2], r[3]) for r in ref])
# (c) mesh refinement
import re as _re
def _lab(c):
    def f(m): return '$a=%g$, $\\lambda=%g$' % (float(m.group(1)), float(m.group(2)))
    c = _re.sub(r'a(\d+\.?\d*)/lam(\d+\.?\d*)', f, c)
    return c.replace(', radius ratio 0.6', ',\nradius ratio 0.6')
M = json.load(open(os.path.join(LIT, 'ripple_mesh_test.json'))); ax = axs[2]
cases = []; [cases.append(m['case']) for m in M if m['case'] not in cases]; rows = []
cm = plt.cm.Greens(np.linspace(0.4, 0.95, len(cases)))
for c, col in zip(cases, cm):
    q = [m for m in M if m['case'] == c and m['h_mesh'] == 0.025]
    if not q: q = [m for m in M if m['case'] == c][-1:]
    q = q[0]; hs = [8, 16, 32]; v = [q[f'prod_h{h}']['ratio_to_mesh'] for h in hs]
    ax.plot([1 / 8, 1 / 16, 1 / 32], v, 'o-', ms=4.5, color=col, lw=1, label=_lab(c))
    for h, x in zip(hs, v): rows.append((c, h, round(x, 5)))
ax.axhline(1, color=COL['base'], lw=1.0); ax.set_xscale('log'); ax.set_xticks([1 / 8, 1 / 16, 1 / 32]); ax.set_xticklabels(['$R/8$', '$R/16$', '$R/32$']); ax.minorticks_off()
ax.set_ylabel('Ratio of contact volume to fine-mesh value'); ax.set_xlabel('production mesh spacing'); ax.set_ylim(0.9895, 1.0008); ax.legend(frameon=False, loc='lower left', ncol=1, handlelength=1.5)
write_csv('figS7c_mesh_refinement', ['case', 'spacing_denominator', 'ratio_to_fine_mesh'], rows)
for lab, a_ in zip('abc', axs): panel_label(a_, lab, dx=-0.25 if lab != 'b' else -0.45)
savefig(fig, 'figS7_checks'); print('done')
