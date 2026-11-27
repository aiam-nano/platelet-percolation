from fig_common import *
from figgen import *; from render import *
L = 16.0; VIEW = (-28, 24); Y0, Y1 = 6.5, 9.2
SCEN = {
 'thin':   dict(kind='thin', A=1e9, seed=11, label='zero-thickness discs'),
 'cyl10':  dict(kind='cyl', A=10, seed=12, label='$A=10$'),
 'cyl100': dict(kind='cyl', A=100, seed=13, label='$A=100$'),
 'stack3': dict(kind='cyl', A=30, N=3, seed=14, label='$A=30$, $N=3$'),
 'stack5': dict(kind='cyl', A=50, N=5, seed=15, label='$A=50$, $N=5$'),
 'ln3':    dict(kind='cyl', A=30, sig=0.3, seed=16, label='$\\sigma=0.3$'),
 'ln5':    dict(kind='cyl', A=30, sig=0.5, seed=17, label='$\\sigma=0.5$'),
 'al_a':   dict(kind='cyl', A=30, ell=4.0, kappa=4, seed=18, label='$d=4$, weak'),
 'al_b':   dict(kind='cyl', A=30, ell=4.0, kappa=20, seed=19, label='$d=4$, strong'),
 'rip1':   dict(kind='ripple', amp_slope=1.0, lam=1.0, seed=20, label='$a=1$, $\\lambda=1$'),
 'rip3':   dict(kind='ripple', amp_slope=3.0, lam=1.0, seed=21, label='$a=3$, $\\lambda=1$'),
}
_cache = {}
def scene(key):
    if key in _cache: return _cache[key]
    sp = dict(SCEN[key]); lab = sp.pop('label'); cfg = make_cfg(sp.pop('kind'), L=L, **sp)
    cnt, l, h, wr = run_trace(cfg)
    info = dict(cfg=cfg, count=cnt, labels=l, hist=h, root=wr, label=lab)
    nrm = cfg['rot'][:cnt, :, 2]; info['S'] = pair_order(cfg['pos'][:cnt], nrm, L) if key.startswith('al') else None
    A = cfg['A']; Rm = cfg['rad'][:cnt]; t = 2 * cfg['hth'][:cnt]
    info['phi'] = cnt / L ** 3 * math.pi * np.mean(Rm ** 2) * np.mean(t) * 100 if cfg['kind'] != 'ripple' else None
    info['span'] = int((l == wr).sum()); _cache[key] = info; return info
def draw_scene(ax, key, nseg=14, title=True, fs=6.8):
    inf = scene(key); cfg = inf['cfg']
    ids = slab_ids(cfg, inf['count'], Y0, Y1)
    allp, r, u, v = draw(ax, cfg, inf['labels'], inf['root'], ids, view=VIEW, nseg=nseg, lw=(0.6 if cfg['kind']=='ripple' else 1.0), edge_scale=(0.9 if cfg['kind']=='ripple' else 0.8))
    b = frame_box(ax, L, Y0, Y1, VIEW)
    ax.set_xlim(b[:, 0].min() - .4, b[:, 0].max() + .4); ax.set_ylim(b[:, 1].min() - .4, b[:, 1].max() + .4); ax.set_aspect('equal'); ax.axis('off')
    if title: ax.set_title(inf['label'], fontsize=fs, loc='center')
    return inf
if __name__ == '__main__':
    import time
    for k in SCEN:
        t = time.time(); i = scene(k); print(k, i['count'], i['span'], i['S'], round(time.time() - t, 1), flush=True)
