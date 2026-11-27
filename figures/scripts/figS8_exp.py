from fig_common import *
import csv
R = list(csv.DictReader(open(os.path.join(LIT, 'experimental_comparison.csv'))))
kc = {'rheological (rigidity)': COL['stack'], 'rheological and electrical': '#56B4E9', 'electrical': COL['base']}
fig, axs = plt.subplots(1, 2, figsize=(7.2, 4.4), gridspec_kw=dict(width_ratios=[1.6, 1], wspace=0.55))
ax = axs[0]; n = len(R); y = np.arange(n)[::-1]; rows = []
ax.axvspan(0.7, 1.3, color=COL['grey'], alpha=0.25, lw=0); ax.axvline(1, color=COL['base'], lw=1.0)
for yi, r in zip(y, R):
    mid, lo, hi = float(r['ratio_mid']), float(r['ratio_min']), float(r['ratio_max']); seg = 'segregated' in r['note']; ex = False
    col = COL['bad'] if ex else kc[r['kind']]
    ax.plot([lo, hi], [yi, yi], color=col, lw=1.1); ax.plot(mid, yi, 's' if seg else 'o', ms=4, color=col, mfc='white' if seg else col, mew=1)
    rows.append((r['system'], r['source'], r['A_low'], r['A_high'], r['sigma_lnR'], mid, lo, hi, r['kind'], 'segregated network' if seg else ('excluded (flake breakage)' if ex else '')))
ax.set_xscale('log'); ax.set_yticks(y); ax.set_yticklabels([r['system'].replace(' (powder coated)', ', powder-coated').replace(' (melt mixed)', ', melt-mixed').split(' (')[0] for r in R]); ax.set_xlabel('Ratio of measured threshold to model'); ax.set_ylim(-0.7, n - 0.3)
ax.set_xlim(0.08, 300); ax.set_xticks([0.1, 0.3, 1, 3, 10, 100]); ax.set_xticklabels(['0.1', '0.3', '1', '3', '10', '100']); ax.minorticks_off()
import matplotlib.lines as ml
fig.legend(handles=[ml.Line2D([], [], marker='o', ls='', color=c, ms=3.8, label=k) for k, c in kc.items()] + [ml.Line2D([], [], marker='s', ls='', mfc='white', color=COL['base'], ms=3.8, label='segregated network')], frameon=False, loc='upper center', bbox_to_anchor=(0.45, 0.03), ncol=4, handletextpad=0.3, columnspacing=1.4)
write_csv('figS8a_experimental_ratios', ['system', 'source', 'A_low', 'A_high', 'sigma', 'ratio_mid', 'ratio_min', 'ratio_max', 'kind', 'flag'], rows)
ax = axs[1]
for r in R:
    ex = False; col = kc[r['kind']]
    lo, hi = float(r['measured_low_vol%']), float(r['measured_high_vol%']); lw_ = float(r['law_vol%_at_A_high']); lw2 = float(r['law_vol%_at_A_low']); lm = math.sqrt(lw_ * lw2); mm = math.sqrt(lo * hi)
    ax.plot([lw_, lw2], [mm, mm], color=col, lw=1.0); ax.plot([lm, lm], [lo, hi], color=col, lw=1.0); ax.plot(lm, mm, 'o', ms=3.4, color=col)
xx = np.logspace(-2.2, 1.3, 10); ax.plot(xx, xx, color=COL['base'], lw=1.0); ax.fill_between(xx, xx * 0.7, xx * 1.3, color=COL['grey'], alpha=0.25, lw=0)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('model threshold (vol. %)'); ax.set_ylabel('measured threshold (vol. %)')
panel_label(axs[0], 'a', dx=-0.6); panel_label(axs[1], 'b', dx=-0.28)
savefig(fig, 'figS8_experimental'); print('ok')
