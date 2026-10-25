"""Open-boundary 'spanning' estimator in the style of Yi and Tawerghi (2009): plates are placed with their centres in a
unit cell with hard (non-periodic) faces, a plate is connected to a face if it crosses the face plane, and the
threshold is the number of plates at which a cluster first joins two opposite faces.  The first-spanning count is
recorded separately for the three axes.  Shapes: zero-thickness disc (mode 0), square or triangle (mode 4).
usage: python3 span_open.py shape L nreal seed"""
import sys, time, math, json
import numpy as np
from numba import njit
from sim2e import disc_contact, poly_contact, random_rotations, dot
from poly_test import polygon

@njit(cache=True)
def uf_find(par, x):
    r = x
    while par[r] != r: r = par[r]
    while par[x] != r:
        nx = par[x]; par[x] = r; x = nx
    return r

@njit(cache=True)
def first_span(mode, pos, rot, L, ncell, br, V, Rdisc, ext):
    """returns int array (3,) : number inserted when axis a first spans (-1 if never); ext = half-extent function:
    for plates the extent along axis a is computed from the vertices."""
    N = pos.shape[0]
    cs = L / ncell
    head = -np.ones((ncell, ncell, ncell), np.int64)
    nxt = -np.ones(N, np.int64)
    par = np.arange(N)
    mask = np.zeros(N, np.int64)     # bit 2a: cluster touches low face of axis a, bit 2a+1: high face
    res = -np.ones(3, np.int64)
    ni = np.empty(3); nj = np.empty(3); zero = np.zeros(3)
    nv = V.shape[0]
    for i in range(N):
        cx = int(pos[i, 0] / cs); cy = int(pos[i, 1] / cs); cz = int(pos[i, 2] / cs)
        cx = min(cx, ncell - 1); cy = min(cy, ncell - 1); cz = min(cz, ncell - 1)
        for q in range(3): ni[q] = rot[i, q, 2]
        # face contacts: extent of the plate along each axis
        for a in range(3):
            lo = 1e300; hi = -1e300
            if mode == 0:
                e = Rdisc * math.sqrt(max(0.0, 1.0 - ni[a] * ni[a])); lo = -e; hi = e
            else:
                for k in range(nv):
                    c = rot[i, a, 0] * V[k, 0] + rot[i, a, 1] * V[k, 1]
                    if c < lo: lo = c
                    if c > hi: hi = c
            if pos[i, a] + lo <= 0.0: mask[i] |= (1 << (2 * a))
            if pos[i, a] + hi >= L: mask[i] |= (1 << (2 * a + 1))
        for dx in range(-1, 2):
            x = cx + dx
            if x < 0 or x >= ncell: continue
            for dy in range(-1, 2):
                y = cy + dy
                if y < 0 or y >= ncell: continue
                for dz in range(-1, 2):
                    z = cz + dz
                    if z < 0 or z >= ncell: continue
                    j = head[x, y, z]
                    while j != -1:
                        dv = pos[j] - pos[i]
                        if dv[0] ** 2 + dv[1] ** 2 + dv[2] ** 2 < (2 * br) ** 2:
                            for q in range(3): nj[q] = rot[j, q, 2]
                            if mode == 0: hit = disc_contact(ni, Rdisc, dv, nj, Rdisc)
                            else: hit = poly_contact(rot[i], dv, rot[j], V)
                            if hit:
                                ra = uf_find(par, i); rb = uf_find(par, j)
                                if ra != rb:
                                    par[ra] = rb; mask[rb] |= mask[ra]
                        j = nxt[j]
        nxt[i] = head[cx, cy, cz]; head[cx, cy, cz] = i
        mk = mask[uf_find(par, i)]
        for a in range(3):
            if res[a] < 0 and (mk >> (2 * a)) & 3 == 3:
                res[a] = i + 1
        if res[0] >= 0 and res[1] >= 0 and res[2] >= 0: break
    return res

def run(shape, L, nreal, seed, Nfac=0.45):
    rng = np.random.default_rng(seed); out = []
    if shape == 'disc': mode, V, br = 0, np.zeros((1, 3)), 1.0
    else: mode, V = 4, polygon(shape); br = float(np.linalg.norm(V, axis=1).max())
    for _ in range(nreal):
        n = int(Nfac * L ** 3); pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng)
        r = first_span(mode, pos, rot, float(L), max(3, int(L // (2 * br))), br, V, 1.0, 0.0)
        out.append(r * 4 / 3 * math.pi / L ** 3)
    return np.array(out)       # (nreal, 3), negative entries mean never spanned

if __name__ == "__main__":
    shape, L, nreal, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    t0 = time.time(); e = run(shape, L, nreal, seed)
    ok = (e > 0).all(axis=1); e = e[ok]
    per_axis = e.mean(axis=0); anyd = e.min(axis=1)
    row = dict(shape=shape, L=L, n=int(ok.sum()), eta_axis_mean=float(per_axis.mean()), eta_axis_err=float(e.mean(axis=1).std(ddof=1) / math.sqrt(len(e))),
               eta_median=float(np.median(e)), eta_any=float(anyd.mean()), eta_any_err=float(anyd.std(ddof=1) / math.sqrt(len(e))), seed=seed, secs=round(time.time() - t0))
    print(json.dumps(row), flush=True)
    open('span_open.jsonl', 'a').write(json.dumps(row) + '\n')
