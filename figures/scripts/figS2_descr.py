from fig_common import *
from matplotlib.patches import Rectangle, Circle
from matplotlib.collections import LineCollection
import fieldgen, align
fig = plt.figure(figsize=(7.2, 7.2))
gs = fig.add_gridspec(3, 12, hspace=0.55, wspace=1.2, height_ratios=[0.8, 0.85, 1.6])
# (a) tactoid cross-sections
ax = fig.add_subplot(gs[0, 0:6]); x0 = 0
for N in (1, 3, 5):
    for k in range(N):
        ax.add_patch(Rectangle((x0, k * 0.34), 5.0, 0.26, fc=COL['stack'], ec='none', alpha=0.95))
    ax.text(x0 + 2.5, -0.15, f'$N={N}$' if N > 1 else '$N=1$', ha='center', va='top', color=COL['base'])
    x0 += 6.4
ax.set_xlim(-0.3, 19.2); ax.set_ylim(-1.0, 2.3); ax.set_aspect('equal'); ax.axis('off')

# (b) ripple profiles
ax = fig.add_subplot(gs[0, 6:12]); x = np.linspace(-1, 1, 400)
for (a_, lam_, c) in ((1, 1, '#7FD1B9'), (3, 1, '#009E73'), (1.5, 2, '#00553E')):
    z = a_ * lam_ / (2 * math.pi) * np.sin(2 * math.pi * x / lam_)
    ax.plot(x, z, color=c, label=f'$a={a_:g}$, $\\lambda={lam_:g}$')
ax.set_aspect('equal'); ax.set_xlabel('position across the sheet ($R$)'); ax.set_ylabel('height ($R$)'); ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1.03, 0.5), ncol=1, handlelength=1.5)
ax.set_ylim(-1.25, 0.75); ax.set_yticks([-0.5, 0, 0.5])
# (c) random wrinkle fields
cfgs = [('R1', 4, 0.4), ('R2', 4, 0.8), ('R3', 4, 1.2), ('R4', 2, 0.8)]
for i, (nm, p, s) in enumerate(cfgs):
    ax = fig.add_subplot(gs[1, 3 * i:3 * i + 3]); g = fieldgen.make_field(1000 * int(nm[1:]) + 3, s, p)
    n = g.shape[0]; ax_ = -fieldgen.GH + fieldgen.HG * np.arange(n); X, Y = np.meshgrid(ax_, ax_); gm = np.where(X * X + Y * Y <= 1.0, g, np.nan)
    im = ax.contourf(X, Y, gm, levels=np.linspace(-0.6, 0.6, 13), cmap='PRGn_r', extend='both'); ax.contour(X, Y, gm, levels=[0], colors='k', linewidths=1.0)
    ax.add_patch(Circle((0, 0), 1, fill=False, ec=COL['base'], lw=1.0)); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(f'{nm}\nrms slope {s:g}', loc='center', linespacing=1.3)
    
# (d) domain maps
Lb = 16.0; n = 6000; rg = np.random.default_rng(5); pos = rg.random((n, 3)) * Lb
nseed = int(round(Lb ** 3 / 4.0 ** 3)); sd = rg.random((nseed, 3)) * Lb; dr = rg.normal(size=(nseed, 3)); dr /= np.linalg.norm(dr, axis=1)[:, None]
gx = np.linspace(0, Lb, 240); GX, GY = np.meshgrid(gx, gx); P3 = np.stack([GX.ravel(), GY.ravel(), np.full(GX.size, 8.0)], 1)
dd = P3[:, None, :] - sd[None]; dd -= Lb * np.round(dd / Lb); dom = np.argmin((dd ** 2).sum(2), 1).reshape(GX.shape)
az = (np.arctan2(dr[:, 1], dr[:, 0]) % math.pi) / math.pi; hue = az[dom]
dd = pos[:, None, :] - sd[None]; dd -= Lb * np.round(dd / Lb); idx = np.argmin((dd ** 2).sum(2), 1); c0_ = dr[idx]
for j, (kap, nm) in enumerate(((4, 'weak scatter: $\\kappa=4$, $S\\approx0.1$'), (20, 'strong scatter: $\\kappa=20$, $S\\approx0.4$'))):
    ax = fig.add_subplot(gs[2, 6 * j:6 * j + 6]); r2 = np.random.default_rng(7)
    ax.imshow(hue, extent=(0, Lb, 0, Lb), origin='lower', cmap='twilight', vmin=0, vmax=1, alpha=0.30, interpolation='nearest')
    ct = align.fisher(kap, n, r2); st = np.sqrt(np.maximum(0, 1 - ct ** 2)); ph = r2.random(n) * 2 * math.pi
    a_ = np.where(np.abs(c0_[:, [0]]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]])); e1 = np.cross(c0_, a_); e1 /= np.linalg.norm(e1, axis=1)[:, None]; e2 = np.cross(c0_, e1)
    nn = ct[:, None] * c0_ + st[:, None] * (np.cos(ph)[:, None] * e1 + np.sin(ph)[:, None] * e2)
    sel = np.abs(pos[:, 2] - 8) < 0.5; p2 = pos[sel]; q = nn[sel]; ang = np.arctan2(q[:, 1], q[:, 0]); ln = np.hypot(q[:, 0], q[:, 1]) * 0.9
    segs = [[(p[0] - 0.5 * l * math.cos(a), p[1] - 0.5 * l * math.sin(a)), (p[0] + 0.5 * l * math.cos(a), p[1] + 0.5 * l * math.sin(a))] for p, a, l in zip(p2, ang, ln)]
    cs = np.abs((q * c0_[sel]).sum(1)); lc = LineCollection(segs, array=cs, cmap='viridis', clim=(0, 1), linewidths=1.2); ax.add_collection(lc)
    ax.set_xlim(0, Lb); ax.set_ylim(0, Lb); ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(nm, loc='center')
    write_csv(f'figS2d_orientation_k{kap}', ['x', 'y', 'nx', 'ny', 'nz'], [(round(a, 3), round(b, 3), round(c, 3), round(d, 3), round(e, 3)) for (a, b, _), (c, d, e) in zip(p2, q)])
for lab, ax in zip('abcd', [fig.axes[0], fig.axes[1], fig.axes[2], fig.axes[6]]): panel_label(ax, lab, dx=-0.08, dy=1.12)
savefig(fig, 'figS2_descriptions'); print('ok')
