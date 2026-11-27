from fig_common import *
from figgen import *; from render import *
L = 16.0; A = 30.0
import os as _o; LW=float(_o.environ.get('LW','0.5')); ES=float(_o.environ.get('ES','0.7'))
cfg = make_cfg('cyl', L=L, seed=3, A=A)
cnt, lab, hist, wr = run_trace(cfg)
view = (-28, 24); y0, y1 = 5.5, 9.5
fracs = [0.25, 0.5, 0.75, 1.0]
fig = plt.figure(figsize=(7.2, 4.9))
gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.62], hspace=0.12, wspace=0.03)
allb = frame_box(plt.figure(figsize=(1, 1)).gca(), L, y0, y1, view); plt.close('all')
fig = plt.figure(figsize=(7.2, 4.9)); gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 0.62], hspace=0.12, wspace=0.03)
n_slab_tot = []
for k, f in enumerate(fracs):
    n = int(round(f * cnt)); n = cnt if f == 1.0 else n
    c, l, h, w = run_trace(cfg, stop=n) if f < 1.0 else (cnt, lab, hist, wr)
    # largest cluster at this stage
    un, ct = np.unique(l, return_counts=True); big = un[np.argmax(ct)] if f < 1.0 else wr
    ids = slab_ids(cfg, n, y0, y1)
    ax = fig.add_subplot(gs[0, k])
    allp, r, u, v = draw(ax, cfg, l, big, ids, view=view, nseg=14, lw=LW, edge_scale=ES)
    b = frame_box(ax, L, y0, y1, view)
    ax.set_xlim(allb[:, 0].min() - 0.4, allb[:, 0].max() + 0.4); ax.set_ylim(allb[:, 1].min() - 0.4, allb[:, 1].max() + 0.4); ax.set_aspect('equal'); ax.axis('off')
    phi = n / L ** 3 * math.pi * 2 / A * 100
    ax.set_title(f'{int(100 * f)}% of the first-wrap count\n{n} flakes · {phi:.2f} vol. %', fontsize=8, loc='center', linespacing=1.3)
    n_slab_tot.append((f, n, len(ids), int(ct.max()) if f < 1.0 else int((lab == wr).sum())))
ax = fig.add_subplot(gs[1, :]); x = np.arange(1, len(hist) + 1); xs = x / cnt
ax.plot(xs, hist / cnt * 100, color=COL['stack'], lw=1.5)
ax.set_xlabel('Fraction of the first-wrap count'); ax.set_ylabel('largest cluster\n(% of flakes at first wrap)')
for f in fracs: ax.axvline(f, color=COL['grey'], lw=1.0, ls=':')
ax.axvline(1.0, color=COL['sim'], lw=1.0); ax.text(1.0, 0.5 * ax.get_ylim()[1] if False else 62, ' first wrapping\n cluster', fontsize=8, ha='left', va='center')
ax.set_xlim(0, 1.18); ax.set_ylim(0, 75)
write_csv('fig2_buildup_curve', ['flakes_inserted', 'largest_cluster'], [(int(a), int(b)) for a, b in zip(x, hist)])
write_csv('fig2_buildup_frames', ['fraction', 'flakes', 'flakes_in_slab_shown', 'largest_cluster'], n_slab_tot)
panel_label(fig.axes[0], 'a', -0.12, 1.2); panel_label(ax, 'b', -0.07, 1.04)
savefig(fig, os.environ.get('OUT','Fig2_buildup')); print(cnt, n_slab_tot)
