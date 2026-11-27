from fig_common import *
from scipy.special import ellipe
def kx(a):
    th = np.linspace(0, np.pi / 2, 4001); return 1 / np.mean(np.sqrt(1 + a * a * np.cos(th) ** 2))
def Kell(a): k = kx(a); return math.pi / (2 * k * ellipe(1 - k * k))
def phi(A, N=1, sig=0.0, S=0.0, d=4.0, a=0.0, lam=1.0, B=0.663, kap=0.625, pf=0.588, gam=0.576, crip=0.647):
    Ae = A / N; Kal = 1 / (1 - pf * S * (1 + gam / (d / 2) ** 2)) if S > 0 else 1.0
    Kr = 1 + (Kell(a) - 1) * math.exp(-crip * a * lam) if a > 0 else 1.0; K = Kal * Kr
    base = c0 * (1 + B / (Ae + 1.5)) / (Ae + c1 + c2 / Ae); dil = (Ae + c1 + c2 / Ae) / (Ae / K + c1 + c2 / Ae)
    return base * dil * Pdisp(sig, Ae, kap)
print('check', phi(30, 2, .3, .2, 4, 1, 1), LAW.phi(30, 2, .3, .2, 4, 1, 1)); assert abs(phi(30, 2, .3, .2, 4, 1, 1) / LAW.phi(30, 2, .3, .2, 4, 1, 1) - 1) < 1e-4
rg = np.random.default_rng(1); M = 300
draw = dict(B=rg.normal(0.663, 0.043, M), kap=rg.normal(0.625, 0.032, M), pf=rg.normal(0.588, 0.011, M), gam=rg.normal(0.576, 0.075, M), crip=rg.normal(0.647, 0.011, M))
def band(fn, xs, keys):
    cen = np.array([fn(x, {}) for x in xs]); sm = np.array([[fn(x, {k: draw[k][i] for k in keys}) for x in xs] for i in range(M)])
    lo, hi = np.percentile(sm, [16, 84], axis=0); return cen, lo, hi
fig, AX = plt.subplots(2, 4, figsize=(7.2, 4.4), gridspec_kw=dict(wspace=0.55, hspace=0.5, height_ratios=[1.3, 1])); axs = AX[0]; rows = []; rel = []
def plot(ax, xs, cen, lo, hi, col, xl, yl, title, key, logx=False, logy=False):
    ax.fill_between(xs, lo, hi, color=col, alpha=0.3, lw=0); ax.plot(xs, cen, color=col); ax.set_xlabel(xl); ax.set_ylabel(yl)
    if logx: ax.set_xscale('log')
    if logy: ax.set_yscale('log')
    rel.append((xs, 100 * (hi - lo) / 2 / cen, col))
    for x, c, l, h in zip(xs, cen, lo, hi): rows.append((key, round(float(x), 4), round(float(c), 5), round(float(l), 5), round(float(h), 5)))
xs = np.logspace(0.7, 3, 25); c, l, h = band(lambda x, o: 100 * phi(x, **o), xs, ['B']); plot(axs[0], xs, c, l, h, COL['base'], '$A$', '$\\phi_c$ (vol. %)', 'thickness', 'thickness', True, True)
xs = np.linspace(0, 0.55, 20); c, l, h = band(lambda x, o: phi(30, sig=x, **o) / phi(30, sig=0), xs, ['kap']); plot(axs[1], xs, c, l, h, COL['disp'], '$\\sigma$', 'Ratio of $\\phi_c$ to that of\nthe monodisperse system', 'dispersity', 'dispersity')
xs = np.linspace(0.02, 0.5, 20); c, l, h = band(lambda x, o: phi(33.333, S=x, d=3, **o) / phi(33.333), xs, ['pf', 'gam']); plot(axs[2], xs, c, l, h, COL['align'], '$S$', 'Ratio of $\\phi_c$ to that of\nthe isotropic system', 'alignment, $d=3$', 'alignment')
xs = np.linspace(0, 3.5, 20); c, l, h = band(lambda x, o: phi(1e6, a=x, lam=1.0, **o) / phi(1e6), xs, ['crip']); plot(axs[3], xs, c, l, h, COL['ripple'], '$a$ ($\\lambda=1$)', 'Ratio of $\\phi_c$ to that of a flat sheet', 'ripple', 'ripple')
for (xs_, w, col), a_, xl in zip(rel, AX[1], ['$A$', '$\\sigma$', '$S$', '$a$ ($\\lambda=1$)']):
    a_.plot(xs_, w, color=col); a_.set_xlabel(xl); a_.set_ylabel('band half-width (%)'); a_.set_ylim(bottom=0)
    if xl == '$A$': a_.set_xscale('log')
for lab, a_ in zip('abcdefgh', list(AX[0]) + list(AX[1])): panel_label(a_, lab, dx=-0.3)
write_csv('figS6_bands', ['panel', 'x', 'central', 'lower_16pct', 'upper_84pct'], rows)
savefig(fig, 'figS6_uncertainty'); print([ (k, r[3], r[4]) for k, r in zip(['a','b','c','d'], [rows[12], rows[25+10], rows[45+10], rows[65+10]])])
