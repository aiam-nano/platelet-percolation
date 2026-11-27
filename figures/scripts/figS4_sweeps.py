from fig_common import *
import new_preds as NP
fig = plt.figure(figsize=(7.2, 5.9)); gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.55)
# (a) phi_c vs A for several N
ax = fig.add_subplot(gs[0, 0]); A = np.logspace(0.9, 3.0, 60); cN = {1: '#b3d4ea', 2: '#7fb5da', 5: '#3f94c8', 10: '#0072B2', 20: '#003d63'}
for N, c in cN.items():
    Av = A[A / N >= 1.5]; ax.plot(Av, [law(a, N=N) for a in Av], color=c, label=f'$N={N}$')
rows = []
for a, y in FL.base: rows.append((a, 1, 100 * y / FS / a))
ax.plot([r[0] for r in rows], [r[2] for r in rows], 'o', ms=3, color=COL['sim'], mfc='white', mew=1.0, label='simulation, $N=1$')
bs = {}
for a, y in FL.base:
    for n in (2, 5, 10, 20):
        if abs(a - 100 / n) < 0.02: bs[(n, 100)] = 100 * y / FS / a
ax.plot([100] * len(bs), list(bs.values()), 's', ms=4.5, color=COL['sim'], mfc=COL['stack'], mew=1.0, label='simulation, $A=100$')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('aspect ratio of one layer, $A$'); ax.set_ylabel('$\\phi_c$ (vol. %)'); ax.legend(frameon=False, ncol=3, loc='lower left', bbox_to_anchor=(-0.05, 1.01), columnspacing=1.0, handlelength=1.4)
write_csv('figS4a_layers_vs_A', ['A', 'N', 'sim_vol_pct_div_1.015'], rows + [(100, n, v) for (n, _), v in bs.items()])
# (b) S-d plane
ax = fig.add_subplot(gs[0, 1]); al = [x for x in FL.al if x[0] < 100]
dv = np.linspace(1.6, 8.8, 60); Sv = np.linspace(0.05, 0.55, 60)
Z = np.array([[LAW.phi(33.333, S=S_, d=d_) / LAW.phi(33.333) for d_ in dv] for S_ in Sv])
cs = ax.contourf(dv, Sv, Z, levels=np.arange(1.0, 2.2, 0.1), cmap='Purples', alpha=0.85); ax.contour(dv, Sv, Z, levels=np.arange(1.0, 2.2, 0.1), colors='white', linewidths=1.0)
cb = fig.colorbar(cs, ax=ax, pad=0.02, fraction=0.05); cb.set_label('model, $\\phi_c$ relative to isotropic'); cb.ax.tick_params(labelsize=6)
ax.scatter([x[1] for x in al], [x[2] for x in al], c=[x[3] for x in al], cmap='Purples', vmin=1.0, vmax=2.1, s=26, edgecolors=COL['sim'], linewidths=1.0)
for i, x in enumerate(al): ax.text(x[1] + 0.12, x[2] + (0.012 if i % 2 else -0.03), f'{x[3]:.2f}', color=COL['sim'], fontsize=7)
ax.set_xlabel('domain diameter $d$ ($R$)'); ax.set_ylabel('local pair order $S$')
write_csv('figS4b_S_d', ['d', 'S', 'sim_ratio', 'se'], [(x[1], x[2], x[3], x[4]) for x in al])
# (c) a-lambda plane
ax = fig.add_subplot(gs[1, 0]); av = np.linspace(0.05, 6.3, 90); lv = np.linspace(0.3, 2.6, 70)
Z = np.array([[LAW.Krip(a_, l_) for a_ in av] for l_ in lv]); cs = ax.contourf(av, lv, Z, levels=np.arange(1.0, 1.9, 0.1), cmap='Greens', alpha=0.85); ax.contour(av, lv, Z, levels=np.arange(1.0, 1.9, 0.1), colors='white', linewidths=1.0)
cb = fig.colorbar(cs, ax=ax, pad=0.02, fraction=0.05); cb.set_label('model $K_{rip}$'); cb.ax.tick_params(labelsize=6)
r = pickle.load(open(os.path.join(SP, 'ripthin2_G.pkl'), 'rb')); fl = r['flat']; m0 = fl.mean(); pts = []
for k, v in r.items():
    if k == 'flat': continue
    pts.append((k[0], k[1], v.mean() / m0))
sc_ = ax.scatter([p[0] for p in pts], [p[1] for p in pts], c=[p[2] for p in pts], cmap='Greens', vmin=1.0, vmax=1.8, s=30, edgecolors=COL['sim'], linewidths=1.0)
for a_, l_, v in pts: ax.text(a_, l_ + (0.13 if int(round(a_ * 2)) % 2 else -0.2), f'{v:.2f}', color=COL['sim'], ha='center', fontsize=7)
ax.axvspan(3.75, 6.4, color=COL['grey'], alpha=0.25, lw=0); ax.text(5.1, 2.5, 'not claimed', ha='center', color=COL['base'])
ax.set_xlabel('peak slope $a$'); ax.set_ylabel('wavelength $\\lambda$ ($R$)'); ax.set_xlim(0, 6.4); ax.set_ylim(0.2, 2.7)
write_csv('figS4c_a_lambda', ['peak_slope', 'wavelength', 'sim_ratio'], [(a, l, round(v, 4)) for a, l, v in pts])
# (d) tactoid distributions
ax = fig.add_subplot(gs[1, 1]); res = {}
for f in ('tact_X.pkl', 'tact_Y.pkl'): res.update(pickle.load(open(os.path.join(SP, f), 'rb')))
rows = []; ids = ['T1', 'T2', 'T3', 'T4', 'T5']
import finallaw2 as L2
for i, A_, k in NP.TACT:
    ks, w = NP.dist(k); Nb = (w * ks).sum(); Nw = (w * ks ** 2).sum() / Nb; f = NP.eig_factor(A_, ks, w)
    pn = 100 * LAW.phi(A_, N=Nb) * f; pw = 100 * LAW.phi(A_, N=Nw) * f; a = res[i][:, 0] * 100 / FS; sim = a.mean(); se = a.std(ddof=1) / math.sqrt(len(a))
    rows.append((i, A_, k, round(Nb, 3), round(Nw, 3), round(pn, 4), round(pw, 4), round(sim, 4), round(se, 4)))
x = np.arange(5); ax.errorbar(x, [r[7] for r in rows], yerr=[r[8] for r in rows], fmt='o', ms=4, color=COL['sim'], label='simulation', zorder=3)
ax.plot(x - 0.12, [r[5] for r in rows], '_', ms=11, mew=2, color=COL['stack'], label='model, number-average $N$')
ax.plot(x + 0.12, [r[6] for r in rows], '_', ms=11, mew=2, color=COL['bad'], label='model, weight-average $N$')
ax.set_xticks(x); ax.set_xticklabels([f'{r[0]}\n$A$={r[1]}' for r in rows]); ax.set_ylabel('$\\phi_c$ (vol. %)'); ax.set_yscale('log'); ax.legend(frameon=False, loc='upper left'); ax.set_ylim(2, 22)
write_csv('figS4d_tactoid_distributions', ['id', 'A', 'distribution', 'Nbar', 'Nw', 'law_number_avg', 'law_weight_avg', 'sim_div_1.015', 'se'], rows)
for lab, a_ in zip('abcd', [fig.axes[0], fig.axes[1], fig.axes[3], fig.axes[5]]): panel_label(a_, lab, dx=-0.17, dy=1.06)
savefig(fig, 'figS4_sweeps'); print(rows)
