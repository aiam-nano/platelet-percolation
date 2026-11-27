from fig_common import *
from matplotlib.ticker import FixedLocator, FixedFormatter
rows = [('Disc', L, n, eta, err, 0.9614) for L, n, eta, err in [(16,600,.9839,.0029),(24,400,.9803,.0022),(32,300,.9765,.0021),(48,150,.9777,.0017),(64,100,.9765,.0017),(96,40,.9751,.0017),(128,16,.9778,.0017)]]
rows += [('Square', L, n, eta, err, 0.8647) for L, n, eta, err in [(16,600,.8802,.0027),(24,400,.8815,.002),(32,300,.8788,.0019),(48,150,.8774,.0016),(64,100,.8808,.0014),(96,40,.8797,.0015),(128,16,.8798,.0018)]]
rows += [('Triangle', L, n, eta, err, 0.7295) for L, n, eta, err in [(16,600,.7487,.0024),(24,400,.7485,.0021),(32,300,.7495,.0016),(48,150,.7491,.0014),(64,100,.7506,.0012),(96,40,.7496,.0012),(128,16,.7523,.0011)]]
sty = {'Disc': (COL['stack'], 'o'), 'Square': (COL['disp'], 's'), 'Triangle': (COL['align'], '^')}
fig = plt.figure(figsize=(7.2, 3.1)); ax = fig.add_axes([0.085, 0.16, 0.68, 0.80])
out = []; ends = []
for s, (c, m) in sty.items():
    r = [x for x in rows if x[0] == s]; L = np.array([x[1] for x in r]); y = np.array([x[3] / x[5] for x in r]); e = np.array([x[4] / x[5] for x in r])
    ax.errorbar(L, y, yerr=e, fmt='-' + m, color=c, ms=4.5, lw=1.0, elinewidth=1.0, capsize=0, mfc=c, mec=c)
    big = [x for x in r if x[1] >= 32]; w = sum(1 / x[4] ** 2 for x in big); mean = sum(x[3] / x[4] ** 2 for x in big) / w
    pl = mean / r[0][5]; ends.append((s, pl, c))
    out += [(s, x[1], x[2], x[3], x[4], x[5], round(x[3] / x[5], 5)) for x in r]
ax.axhline(1, color=COL['base'], lw=1.0, ls=(0, (4, 3))); ax.text(17, 1.0015, 'published value', va='bottom', ha='left', color=COL['base'])
ax.set_xscale('log', base=2); ax.set_xlim(14, 140); ax.set_ylim(0.995, 1.04)
ax.xaxis.set_major_locator(FixedLocator([16, 24, 32, 48, 64, 96, 128])); ax.xaxis.set_major_formatter(FixedFormatter(['16', '24', '32', '48', '64', '96', '128'])); ax.minorticks_off()
ax.set_xlabel('box side $L$ (equivalent radii)'); ax.set_ylabel('Ratio of simulated to published reduced density')
ax.set_yticks([1.00, 1.01, 1.02, 1.03, 1.04]); ax.set_yticklabels(['1.00', '1.01', '1.02', '1.03', '1.04'])
ends.sort(key=lambda t: t[1], reverse=True); ly = [e[1] for e in ends]
for i in range(1, len(ly)):
    if ly[i - 1] - ly[i] < 0.0042: ly[i] = ly[i - 1] - 0.0042
for (s, v, c), yy in zip(ends, ly):
    ax.annotate(f'{s} +{(v - 1) * 100:.1f} %', xy=(128, v), xytext=(150, yy), textcoords='data', color=c, va='center', ha='left', annotation_clip=False, arrowprops=dict(arrowstyle='-', color=COL['grey'], lw=1.0, shrinkA=2, shrinkB=4))
write_csv('fig3_plate_validation', ['shape', 'box_side_L', 'realisations', 'reduced_density', 'standard_error', 'published', 'ratio'], out)
print([(s, round((v - 1) * 100, 2)) for s, v, c in ends])
savefig(fig, 'Fig3_plates')
