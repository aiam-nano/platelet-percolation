from fig_common import *
R = [('Few-layer graphene in NMP (Barwich 2015)',2.01,1.96,2.06,'rigidity'),('xGn-25 in PLA (Sabzi 2013)',101.13,91.7,111.53,'rheological and electrical'),('N02 in PLA (Sabzi 2013)',24.43,21.5,27.75,'rheological and electrical'),
('GNP in polyester (Lu 2006)',1.06,1.06,1.06,'electrical'),('Foliated graphite in nylon 6 (Weng 2005)',1.22,1.22,1.22,'electrical'),('Graphite sheets in PMMA (Zheng 2002)',1.29,1.15,1.43,'electrical'),
('Graphite nanosheets in PMMA (Chen 2003)',0.49,0.22,1.08,'electrical'),('Graphene nanosheets in TPU (Jan 2016)',4.86,4.86,4.86,'electrical'),('GNP-L in PLA (Gao 2017)',6.24,3.89,10.04,'electrical'),
('GNP-S in PLA (Gao 2017)',3.94,2.34,6.67,'electrical'),('xGnP-1 in PP, powder-coated (Kalaitzidou 2010)',0.14,0.14,0.14,'electrical'),('xGnP-1 in PP, melt-mixed (Kalaitzidou 2010)',11.96,11.28,12.68,'electrical'),
('xGnP-15 in PP, powder-coated (Kalaitzidou 2010)',6.25,6.25,6.25,'electrical'),('xGnP-15 in PP, melt-mixed (Kalaitzidou 2010)',176.71,166.6,187.43,'electrical'),('GNP in PLA, AFM sizes (Lashgari 2025)',0.73,0.63,0.85,'rheological and electrical'),('xGnP-M-15 in epoxy (King 2012)',14.73,14.73,14.73,'electrical')]
R.sort(key=lambda r: -r[1]); n = len(R); STEP = 0.24
Ht = n * STEP + 0.75; fig = plt.figure(figsize=(7.2, Ht)); ax = fig.add_axes([0.40, 0.55 / Ht, 0.57, n * STEP / Ht])
ax.set_xscale('log'); ax.set_xlim(0.1, 300); ax.set_ylim(n - 0.4, -0.6)
ax.axvspan(0.7, 1.3, color='#e6e6e6', lw=0, zorder=0)
for t in (0.1, 0.3, 3, 10, 100): ax.axvline(t, color='#d9d9d9', lw=1.0, zorder=0)
ax.axvline(1, color=COL['base'], lw=1.0, ls=(0, (4, 3)), zorder=1)
for i, (nm, mid, lo, hi, kind) in enumerate(R):
    el = kind == 'electrical'; c = '#6e6e6e' if el else COL['stack']
    if hi > lo: ax.plot([lo, hi], [i, i], color=c, lw=1.25, zorder=2, solid_capstyle='butt')
    ax.plot(mid, i, 'o', ms=5, mfc=c if not el else 'white', mec=c, mew=1.25, zorder=3)
    ax.text(0.1 * 0.62, i, nm, ha='right', va='center', clip_on=False)
ax.set_yticks([]); ax.spines['left'].set_visible(False)
ax.set_xticks([0.1, 0.3, 1, 3, 10, 100]); ax.set_xticklabels(['0.1', '0.3', '1', '3', '10', '100']); ax.minorticks_off()
ax.set_xlabel('Ratio of measured threshold to model (logarithmic scale)')
ax.text(1, -0.95, 'model = 1', ha='center', va='bottom')
ax.plot([], [], 'o', ms=5, mfc=COL['stack'], mec=COL['stack'], label='rigidity or rheological'); ax.plot([], [], 'o', ms=5, mfc='white', mec='#6e6e6e', mew=1.25, label='electrical')
ax.legend(frameon=False, loc='lower right', bbox_to_anchor=(1.0, 0.02), handletextpad=0.3, borderaxespad=0.2)
near = sum(1 for r in R if 0.7 <= r[1] <= 1.3); print(near, len(R))
write_csv('fig6_measured_over_model', ['system', 'ratio_mid', 'ratio_low', 'ratio_high', 'kind'], R)
savefig(fig, 'Fig6_measured'); print(Ht)
