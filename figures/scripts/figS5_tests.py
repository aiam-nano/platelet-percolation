from fig_common import *
import pickle
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw=dict(wspace=0.32, width_ratios=[1, 1.15]))
# (a) sigma x N : model curves, A = 100 layers of thickness; relative to the monodisperse single-layer-equivalent
ax = axs[0]; sg = np.linspace(0, 0.9, 50); cc = {1: '#b3d4ea', 3: '#3f94c8', 10: '#003d63'}
for N, c in cc.items():
    ax.plot(sg, [law(100, N=N, sig=s) / law(100, N=N) for s in sg], color=c, label=f'$N={N}$ ($A_e={100 // N}$)')
ax.plot(sg, [law(100, sig=s) / law(100) for s in sg], color=COL['grey'], ls=':', lw=1, label='single-effect curve for $N=1$')
ax.axvspan(0.55, 0.95, color=COL['grey'], alpha=0.18, lw=0)
# sim point: Set 1 (A=100, N=3, sigma=0.3): sim 3.391 vol%
p = 3.391 / FS; ax.plot([0.3], [p / law(100, N=3)], 's', ms=4, color=COL['sim'], label='simulation, set 1')
ax.set_xlim(0, 0.9); ax.text(0.73, 0.97, 'not claimed', rotation=90, va='top', ha='center', color=COL['base'], transform=ax.get_xaxis_transform()); ax.set_xlabel('$\\sigma$, standard deviation of ln $R$'); ax.set_ylabel('Ratio of $\\phi_c$ to that of\nthe monodisperse system'); ax.legend(frameon=False, loc='lower left')
write_csv('figS5a_sigma_N', ['sigma', 'N1', 'N3', 'N10'], [(round(s, 3),) + tuple(round(law(100, N=n, sig=s) / law(100, N=n), 5) for n in (1, 3, 10)) for s in sg])
# (b) sigma x ripple: P1-P4 ratio of rippled to flat at same dispersity
ax = axs[1]; rip = pickle.load(open(os.path.join(SP, 'gaprip_r.pkl'), 'rb')); cases = {'P1': (0.3, 2.0, 1.0), 'P2': (0.5, 2.0, 1.0), 'P3': (0.4, 1.5, 0.7), 'P4': (0.5, 3.0, 1.5)}; rows = []
for i, (k, (s, a, l)) in enumerate(cases.items()):
    o = rip[k]; r = o[:, 1] / o[:, 0]; m = o[:, 1].mean() / o[:, 0].mean(); se = r.std(ddof=1) / math.sqrt(len(r)); lw_ = LAW.Krip(a, l)
    ax.errorbar(i, m, yerr=se, fmt='o', ms=4, color=COL['sim'], zorder=3, label='simulation' if i == 0 else None); ax.plot(i, lw_, '_', ms=14, mew=2, color=COL['ripple'], label='model, $K_{rip}$' if i == 0 else None)
    rows.append((k, s, a, l, round(lw_, 4), round(m, 4), round(se, 4)))
ax.set_xticks(range(4)); ax.set_xticklabels([f'{k}\n$\\sigma$={v[0]:g}\n$a$={v[1]:g}\n$\\lambda$={v[2]:g}' for k, v in cases.items()]); ax.set_ylabel('Ratio of $\\phi_c$, rippled to flat, same $\\sigma$'); ax.set_ylim(0.95, 1.5)
ax.legend(frameon=False, loc='upper left')
write_csv('figS5b_sigma_ripple', ['case', 'sigma', 'a', 'lambda', 'model', 'sim_ratio_of_means', 'se'], rows)
panel_label(axs[0], 'a', dx=-0.25); panel_label(axs[1], 'b', dx=-0.25); savefig(fig, 'figS5_interactions'); print('ok')
