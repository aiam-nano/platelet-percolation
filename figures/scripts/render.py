"""Painter's-algorithm renderer for simulator configurations (orthographic, shaded polygons, depth-sorted)."""
import math, numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgb

TEAL = np.array(to_rgb('#0072B2')); GREY = np.array(to_rgb('#a6a6a6')); PALE = np.array(to_rgb('#dedede'))
LIGHT = np.array([0.35, -0.55, 0.76]); LIGHT /= np.linalg.norm(LIGHT)

def camera(az_deg, el_deg):
    az, el = math.radians(az_deg), math.radians(el_deg)
    v = np.array([math.cos(el) * math.sin(az), -math.cos(el) * math.cos(az), math.sin(el)])  # towards viewer
    r = np.cross([0, 0, 1.0], v); r /= np.linalg.norm(r)
    u = np.cross(v, r)
    return r, u, v

def _circle(m, R=1.0):
    t = np.linspace(0, 2 * math.pi, m, endpoint=False)
    return np.column_stack([R * np.cos(t), R * np.sin(t)])

def particle_polys(cfg, i, nseg=18):
    """Return list of (vertices in world coords (k,3), normal (3,)) for particle i."""
    R = cfg['rad'][i]; M = cfg['rot'][i]; p = cfg['pos'][i]; kind = cfg['kind']
    ex, ey, ez = M[:, 0], M[:, 1], M[:, 2]
    polys = []
    if kind == 'thin':
        c = _circle(nseg + 6, R)
        polys.append((p + c[:, [0]] * ex + c[:, [1]] * ey, ez))
    elif kind == 'cyl':
        c = _circle(nseg, R); N = int(cfg['N']); h = cfg['hth'][i]
        t = 2 * h / N                         # total thickness / layers
        if N == 1:
            zs = [(-h, h)]
        else:
            lt = 0.76 * t
            zs = [(-h + k * t + 0.12 * t, -h + k * t + 0.12 * t + lt) for k in range(N)]
        for (z0, z1) in zs:
            top = p + c[:, [0]] * ex + c[:, [1]] * ey + z1 * ez
            bot = p + c[:, [0]] * ex + c[:, [1]] * ey + z0 * ez
            polys.append((top, ez)); polys.append((bot, -ez))
            for k in range(nseg):
                k2 = (k + 1) % nseg
                q = np.array([bot[k], bot[k2], top[k2], top[k]])
                mid = (c[k] + c[k2]) / 2
                nn = mid[0] * ex + mid[1] * ey
                nn = nn / (np.linalg.norm(nn) + 1e-12)
                polys.append((q, nn))
    else:  # rippled sheet, polar mesh
        sh = cfg['sheet']; amp = cfg['amp']; kk = cfg['kk']
        rings = [0.0, 0.34, 0.67, 1.0]; ns = 20
        def pt(r, th):
            u = r * math.cos(th); v = r * math.sin(th)
            x = np.interp(u, sh['utab'], sh['xtab']); z = amp * math.sin(kk * x)
            return p + x * ex + v * ey + z * ez
        for a in range(len(rings) - 1):
            for k in range(ns):
                th0, th1 = 2 * math.pi * k / ns, 2 * math.pi * (k + 1) / ns
                if a == 0:
                    q = np.array([pt(0, 0), pt(rings[1], th0), pt(rings[1], th1)])
                else:
                    q = np.array([pt(rings[a], th0), pt(rings[a + 1], th0), pt(rings[a + 1], th1), pt(rings[a], th1)])
                n = np.cross(q[1] - q[0], q[2] - q[0]); n /= np.linalg.norm(n) + 1e-12
                polys.append((q, n))
    return polys

def draw(ax, cfg, labels, root_hl, ids, view=(-28, 24), nseg=16, edge=0.12, lw=0.25, palette=None, zoom_pad=0.4, hl_extra=None, edge_scale=0.55):
    """Draw particles `ids` (indices into cfg); particles whose label equals root_hl are teal; others grey."""
    r, u, v = camera(*view)
    P = []; C = []; D = []
    sizes = {}
    if labels is not None:
        un, cnt = np.unique(labels, return_counts=True); sizes = dict(zip(un.tolist(), cnt.tolist()))
    for i in ids:
        lab = labels[i] if labels is not None else -2
        if lab == root_hl: base = TEAL
        elif sizes.get(lab, 1) > 1: base = GREY
        else: base = PALE
        if palette is not None and hl_extra is not None and i in hl_extra: base = palette
        for verts, n in particle_polys(cfg, i, nseg):
            if np.dot(n, v) < 0 and cfg['kind'] == 'cyl':
                continue   # back-facing cap/side of a closed solid
            shade = 0.58 + 0.42 * abs(np.dot(n, LIGHT))
            col = np.clip(base * shade + (1 - shade) * 0.02, 0, 1)
            sx = verts @ r; sy = verts @ u
            P.append(np.column_stack([sx, sy])); C.append(col); D.append((verts.mean(0) @ v))
    order = np.argsort(D)
    P = [P[k] for k in order]; C = [C[k] for k in order]
    ec = [np.clip(c * edge_scale, 0, 1) for c in C]
    pc = PolyCollection(P, facecolors=C, edgecolors=ec, linewidths=lw)
    ax.add_collection(pc)
    allp = np.vstack(P) if P else np.zeros((1, 2))
    return allp, r, u, v

def frame_box(ax, L, y0, y1, view, color='#7f7f7f', lw=1.0):
    """Outline of the slab 0<x<L, y0<y<y1, 0<z<L."""
    r, u, v = camera(*view)
    xs = [0, L]; ys = [y0, y1]; zs = [0, L]
    pts = {(a, b, c): np.array([xs[a], ys[b], zs[c]]) for a in (0, 1) for b in (0, 1) for c in (0, 1)}
    edges = []
    for a in (0, 1):
        for b in (0, 1):
            edges.append(((a, b, 0), (a, b, 1))); edges.append(((a, 0, b), (a, 1, b))); edges.append(((0, a, b), (1, a, b)))
    allp = []
    for e0, e1 in edges:
        p0, p1 = pts[e0], pts[e1]
        ax.plot([p0 @ r, p1 @ r], [p0 @ u, p1 @ u], color=color, lw=lw, zorder=0, solid_capstyle='round')
        allp += [[p0 @ r, p0 @ u], [p1 @ r, p1 @ u]]
    return np.array(allp)

def slab_ids(cfg, n_in, y0, y1):
    p = cfg['pos'][:n_in]
    return np.where((p[:, 1] >= y0) & (p[:, 1] < y1))[0]
