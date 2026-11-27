from fig_common import *
import csv
rows = list(csv.DictReader(open(os.path.join(LIT, 'experimental_comparison.csv'))))
fig = plt.figure(figsize=(7.2, 3.7))
ax = fig.add_axes([0.085, 0.15, 0.43, 0.68]); A = np.logspace(1.0, 4.2, 120)
ideal = lambda x: 100 * c0 / x            # thin-disc threshold 1.44/A in vol. %
for N in (1, 2, 5, 10):
    y = [law(x, N=N) / ideal(x) if x / N > 3 else np.nan for x in A]
    ax.plot(A, y, color=COL['stack'], lw=1.25)
    ax.text(2.6e4, N, f'$N={N}$', color=COL['stack'], va='center')
ax.hlines(1, 14, 1.65e4, color=COL['base'], ls=(0, (4, 3)), lw=1.0)
band_hi = 1.86 * 1.79
ax.fill_between(A, [law(x) / ideal(x) for x in A], [law(x) / ideal(x) * band_hi for x in A], color='#e6e6e6', lw=0, zorder=0)
out = []
for r in rows:
    Al, Ah = float(r['A_low']), float(r['A_high']); lo, hi = float(r['measured_low_vol%']), float(r['measured_high_vol%'])
    am = math.sqrt(Al * Ah); mm = math.sqrt(lo * hi)
    if am > 12000: continue
    rh = 'rheolog' in r['kind']; f = 1 / ideal(am)
    ax.errorbar(am, mm * f, xerr=[[am - Al], [Ah - am]], yerr=[[(mm - lo) * f], [(hi - mm) * f]], fmt='o', ms=4.5, color='#000000' if rh else '#555555', mfc='#000000' if rh else 'white', mew=1.0, lw=1.0, zorder=5)
    out.append((r['system'], round(am, 1), round(mm, 4), round(mm * f, 3), r['kind']))
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(14, 3.5e4); ax.set_ylim(0.04, 400)
ax.set_yticks([0.1, 0.3, 1, 3, 10, 30, 100, 300]); ax.set_yticklabels(['0.1', '0.3', '1', '3', '10', '30', '100', '300']); ax.minorticks_off()
ax.set_xlabel('aspect ratio of the flakes $A=2R/t$'); ax.set_ylabel('Ratio of measured threshold to\nideal thin-disc value ($1.44/A$)')
ax.plot([], [], 'o', ms=4.5, color='#000000', label='measured, rheological'); ax.plot([], [], 'o', ms=4.5, color='#555555', mfc='white', label='measured, electrical')
ax.plot([], [], color=COL['stack'], label='model, $N$ layers per flake'); ax.fill_between([], [], [], color='#e6e6e6', label='model for single layers with the most\nalignment and ripple tested'); ax.legend(frameon=False, loc='lower left', bbox_to_anchor=(-0.17, 1.02), ncol=2, fontsize=8, handlelength=1.5, columnspacing=1.2, labelspacing=0.3)
panel_label(ax, 'a', -0.17, 1.30)
write_csv('fig4a_measured', ['system', 'A_geomean', 'measured_geomean_vol_pct', 'measured_over_1.44_per_A', 'kind'], out)
# (b) size of each lever over the range tested ------------------------------------------------
ax = fig.add_axes([0.765, 0.15, 0.20, 0.68]); A0 = 100.0
lev = [('restacking\n$N$ = 1 to 10', None, law(A0, N=10) / law(A0)),
       ('alignment\n$S$ = 0.5, $d$ = 2', None, LAW.phi(A0, S=0.5, d=2) / LAW.phi(A0)),
       ('ripple\n$a$ = 3.5, $\\lambda$ = 0.5', None, LAW.phi(A0, a=3.5, lam=0.5) / LAW.phi(A0)),
       ('thickness\n$A$ = 10 vs thin limit', None, law(10) / ideal(10)),
       ('size dispersity\n$\\sigma$ = 0.5', None, law(A0, sig=0.5) / law(A0))]
ys = np.arange(len(lev))[::-1]
for y, (nm, c, v) in zip(ys, lev):
    ax.barh(y, math.log10(v), left=0, color=COL['stack'] if v > 1 else COL['disp'], height=0.55)
    ax.text(1.05, y, f'×{v:.2f}', va='center', ha='left', transform=ax.get_yaxis_transform())
ax.axvline(0, color=COL['sim'], lw=1.0); ax.set_yticks(ys); ax.set_yticklabels([x[0] for x in lev], fontsize=8, linespacing=1.15)
ax.set_xlim(-0.3, 1.0); ax.set_xticks(np.log10([0.5, 1, 2, 5, 10])); ax.set_xticklabels(['0.5', '1', '2', '5', '10'])
ax.set_xlabel('Ratio of $\\phi_c$ with to without\nthe departure (log scale)'); panel_label(ax, 'b', -1.25, 1.30)
write_csv('fig4b_levers', ['lever', 'factor_on_phi_c'], [(x[0].replace('\n', ' '), round(float(x[2]), 4)) for x in lev])
savefig(fig, 'Fig7_design_map'); print([(x[0].replace('\n', ' '), round(float(x[2]), 3)) for x in lev])
