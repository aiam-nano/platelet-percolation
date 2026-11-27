from fig_common import *
import align
Lb = 24.0; ell = 6.0
def make(kappa, seed, n=7000):
    rg = np.random.default_rng(seed); pos = rg.random((n, 3)) * Lb
    nseed = int(round(Lb ** 3 / ell ** 3)); sd = rg.random((nseed, 3)) * Lb; dr = rg.normal(size=(nseed, 3)); dr /= np.linalg.norm(dr, axis=1)[:, None]
    d = pos[:, None, :] - sd[None]; d -= Lb * np.round(d / Lb); idx = np.argmin((d ** 2).sum(2), 1); c0_ = dr[idx]
    ct = align.fisher(kappa, n, rg) if kappa > 0 else 1 - 2 * rg.random(n); st = np.sqrt(np.maximum(0, 1 - ct ** 2)); ph = rg.random(n) * 2 * math.pi
    a_ = np.where(np.abs(c0_[:, [0]]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]])); e1 = np.cross(c0_, a_); e1 /= np.linalg.norm(e1, axis=1)[:, None]; e2 = np.cross(c0_, e1)
    nn = ct[:, None] * c0_ + st[:, None] * (np.cos(ph)[:, None] * e1 + np.sin(ph)[:, None] * e2)
    Sd = float(np.mean(1.5 * ct ** 2 - 0.5)); Sp = align.shell_P2(pos, nn, Lb, 2.0)[0]
    return ct, Sd, Sp
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw=dict(wspace=0.32))
ax = axs[0]; cc = {4: '#CC79A7', 20: '#5E2146'}; x = np.linspace(-1, 1, 400)
for kap in (4, 20):
    ct, Sd, Sp = make(kap, 3, 20000)
    ax.hist(ct, bins=40, range=(-1, 1), density=True, color=cc[kap], alpha=0.35, lw=0)
    f = kap * np.exp(kap * x) / (np.exp(kap) - np.exp(-kap)) * 0 + kap * np.exp(kap * (x - 1)) / (1 - np.exp(-2 * kap))
    ax.plot(x, f, color=cc[kap], label=f'$\\kappa={kap}$ (Fisher density)')
ax.set_xlabel('$\\cos\\theta$ between normal and domain director'); ax.set_ylabel('probability density'); ax.legend(frameon=False, loc='upper left')
ax.set_yscale('log'); ax.set_ylim(0.05, 40)
ax = axs[1]; ks = [0, 1, 2, 3, 4.733, 6.312, 10, 20, 40]; rows = []
for kap in ks:
    for sd_ in (1, 2, 3):
        ct, Sd, Sp = make(kap, int(10 * sd_ + kap), 7000); rows.append((kap, sd_, Sd, Sp))
Sd = np.array([r[2] for r in rows]); Sp = np.array([r[3] for r in rows])
ax.plot(Sd, Sp, 'o', ms=4.5, color=COL['align'], mfc='white', mew=1.0, label='simulated generator')
m = float((Sd ** 2 * Sp).sum() / (Sd ** 4).sum()); xx = np.linspace(0, 1, 60); ax.plot(xx, m * xx ** 2, color=COL['align'], label=f'fit: {m:.3f} $\\times$ (director order)$^2$')
ax.plot(xx, 0.665 * xx ** 2, color=COL['grey'], ls='--', lw=1, label='0.665 used in the model')
ax.plot(xx, xx, color=COL['grey'], ls=':', lw=1, label='one-to-one')
ax.set_xlabel('director order $\\langle P_2(\\cos\\theta)\\rangle$'); ax.set_ylabel('pair order $S$ (neighbours within $2R$)'); ax.legend(frameon=False, loc='upper left'); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
write_csv('figS3b_orderMapping', ['kappa', 'seed', 'director_order', 'pair_order'], [(a, b, round(c, 4), round(d, 4)) for a, b, c, d in rows])
for lab, a_ in zip('ab', axs): panel_label(a_, lab, dx=-0.12)
savefig(fig, 'figS3_orientation'); print('coef', m, [(r[0], round(r[2],3), round(r[3],3)) for r in rows if r[1]==1])
