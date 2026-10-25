"""Independent test of the curved-sheet (ripple) contact used for the frozen results.

Reference method: each sheet is represented by a fine triangular mesh (Delaunay triangulation of a lattice of
material points plus the rim, mapped onto z = A sin(kx) with the arc-length coordinate obtained from the
closed-form elliptic integral), and two sheets are in contact when any triangle of one intersects any triangle of
the other (exact segment-triangle tests).  Nothing is shared with the production test (sim4.sheet_contact), which
uses a coarse edge mesh of one sheet crossed against the implicit surface of the other, interpolated arc-length
table and an analytic domain test.

usage: python3 ripple_mesh_test.py npairs seed"""
import sys, math, time, json
import numpy as np
from numba import njit
from scipy.spatial import Delaunay
from scipy.special import ellipeinc
import sim4

# ---------------------------------------------------------------- mesh construction (independent of sim4.make_sheet)
def arclength_table(amp, lam, xmax, n=400001):
    k = 2 * math.pi / lam; B = amp * k
    x = np.linspace(-xmax, xmax, n); t = k * x
    # closed form: int_0^t sqrt(1+B^2 cos^2 s) ds = E(pi/2 | -B^2) - E(pi/2 - t | -B^2)
    F = ellipeinc(math.pi / 2, -B * B) - ellipeinc(math.pi / 2 - t, -B * B)
    return x, F / k

def fine_mesh(amp, lam, h, R0=1.0):
    """triangles (nt,3,3) of the rippled sheet whose material disc has radius R0; lattice spacing h."""
    pts = []
    # hexagonal lattice in material coordinates (u,v)
    nj = int(R0 / (h * math.sqrt(3) / 2)) + 2
    for j in range(-nj, nj + 1):
        for i in range(-int(R0 / h) - 2, int(R0 / h) + 3):
            u = i * h + (0.5 * h if j % 2 else 0.0); v = j * h * math.sqrt(3) / 2
            if u * u + v * v < (R0 - 0.5 * h) ** 2: pts.append((u, v))
    nrim = int(math.ceil(2 * math.pi * R0 / h))
    for k in range(nrim):
        a = 2 * math.pi * k / nrim; pts.append((R0 * math.cos(a), R0 * math.sin(a)))
    P = np.array(pts); tri = Delaunay(P).simplices
    if amp == 0:
        X = np.column_stack([P[:, 0], P[:, 1], np.zeros(len(P))])
    else:
        xs, us = arclength_table(amp, lam, 1.6 * R0)
        x = np.interp(P[:, 0], us, xs)
        X = np.column_stack([x, P[:, 1], amp * np.sin(2 * math.pi / lam * x)])
    T = X[tri]
    # discard sliver triangles outside the disc (Delaunay of a convex polygon cannot create them), keep all
    return T, X

# ---------------------------------------------------------------- geometry
@njit(cache=True)
def seg_tri(p, q, a, b, c):
    """segment pq crosses triangle abc (Moller-Trumbore, closed)"""
    e1 = b - a; e2 = c - a; d = q - p
    h = np.cross(d, e2); det = e1[0] * h[0] + e1[1] * h[1] + e1[2] * h[2]
    if abs(det) < 1e-14: return False
    f = 1.0 / det; s = p - a
    u = f * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2])
    if u < 0.0 or u > 1.0: return False
    qv = np.cross(s, e1)
    v = f * (d[0] * qv[0] + d[1] * qv[1] + d[2] * qv[2])
    if v < 0.0 or u + v > 1.0: return False
    t = f * (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2])
    return 0.0 <= t <= 1.0

@njit(cache=True)
def tri_tri(A, B):
    for k in range(3):
        if seg_tri(A[k], A[(k + 1) % 3], B[0], B[1], B[2]): return True
    for k in range(3):
        if seg_tri(B[k], B[(k + 1) % 3], A[0], A[1], A[2]): return True
    return False

@njit(cache=True)
def mesh_contact(TA, TB, cs, lo, nc):
    """any triangle of TA meets any triangle of TB; TB binned in a uniform grid."""
    nb = TB.shape[0]
    cnt = np.zeros(nc * nc * nc + 1, np.int64)
    bmin = np.empty((nb, 3)); bmax = np.empty((nb, 3))
    for t in range(nb):
        for q in range(3):
            bmin[t, q] = min(TB[t, 0, q], TB[t, 1, q], TB[t, 2, q]); bmax[t, q] = max(TB[t, 0, q], TB[t, 1, q], TB[t, 2, q])
    for t in range(nb):
        i0 = max(0, int((bmin[t, 0] - lo) / cs)); i1 = min(nc - 1, int((bmax[t, 0] - lo) / cs))
        j0 = max(0, int((bmin[t, 1] - lo) / cs)); j1 = min(nc - 1, int((bmax[t, 1] - lo) / cs))
        k0 = max(0, int((bmin[t, 2] - lo) / cs)); k1 = min(nc - 1, int((bmax[t, 2] - lo) / cs))
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                for k in range(k0, k1 + 1):
                    cnt[(i * nc + j) * nc + k + 1] += 1
    for c in range(nc * nc * nc): cnt[c + 1] += cnt[c]
    fill = cnt[:-1].copy(); items = np.empty(cnt[-1], np.int64)
    for t in range(nb):
        i0 = max(0, int((bmin[t, 0] - lo) / cs)); i1 = min(nc - 1, int((bmax[t, 0] - lo) / cs))
        j0 = max(0, int((bmin[t, 1] - lo) / cs)); j1 = min(nc - 1, int((bmax[t, 1] - lo) / cs))
        k0 = max(0, int((bmin[t, 2] - lo) / cs)); k1 = min(nc - 1, int((bmax[t, 2] - lo) / cs))
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                for k in range(k0, k1 + 1):
                    c = (i * nc + j) * nc + k; items[fill[c]] = t; fill[c] += 1
    for a in range(TA.shape[0]):
        amin0 = min(TA[a, 0, 0], TA[a, 1, 0], TA[a, 2, 0]); amax0 = max(TA[a, 0, 0], TA[a, 1, 0], TA[a, 2, 0])
        amin1 = min(TA[a, 0, 1], TA[a, 1, 1], TA[a, 2, 1]); amax1 = max(TA[a, 0, 1], TA[a, 1, 1], TA[a, 2, 1])
        amin2 = min(TA[a, 0, 2], TA[a, 1, 2], TA[a, 2, 2]); amax2 = max(TA[a, 0, 2], TA[a, 1, 2], TA[a, 2, 2])
        i0 = max(0, int((amin0 - lo) / cs)); i1 = min(nc - 1, int((amax0 - lo) / cs))
        j0 = max(0, int((amin1 - lo) / cs)); j1 = min(nc - 1, int((amax1 - lo) / cs))
        k0 = max(0, int((amin2 - lo) / cs)); k1 = min(nc - 1, int((amax2 - lo) / cs))
        if i1 < i0 or j1 < j0 or k1 < k0: continue
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                for k in range(k0, k1 + 1):
                    c = (i * nc + j) * nc + k
                    for p in range(cnt[c], cnt[c + 1]):
                        t = items[p]
                        if bmin[t, 0] > amax0 or bmax[t, 0] < amin0 or bmin[t, 1] > amax1 or bmax[t, 1] < amin1 or bmin[t, 2] > amax2 or bmax[t, 2] < amin2: continue
                        if tri_tri(TA[a], TB[t]): return True
    return False

@njit(cache=True)
def rotate_tris(T, R, shift, out):
    for t in range(T.shape[0]):
        for v in range(3):
            for q in range(3):
                out[t, v, q] = R[q, 0] * T[t, v, 0] + R[q, 1] * T[t, v, 1] + R[q, 2] * T[t, v, 2] + shift[q]

def sim_rand_rot(n, rng):
    return sim4.random_rotations(n, rng)

def run(amp, lam, h_mesh, npairs, seed, sr=1.0, prod_h=(8, 16, 32)):
    """sr = radius of sheet j relative to sheet i (polydisperse pair)"""
    rng = np.random.default_rng(seed)
    kind = 'ripple' if amp > 0 else 'flat'
    shs = {q: sim4.make_sheet(kind, 1.0, h=1 / q, amp=amp, lam=lam) for q in prod_h}   # production meshes: h = R/8 is the one used for all frozen results
    T, X = fine_mesh(amp, lam, h_mesh)
    brad = float(np.linalg.norm(X, axis=1).max()); rmax = brad * (1 + sr)
    offs = np.zeros(1)
    TAr = np.empty_like(T); TBr = np.empty_like(T)
    cs = 4 * h_mesh + 0.05; lo = -(rmax + 0.2); nc = int(2 * (rmax + 0.2) / cs) + 1
    res = []
    for _ in range(npairs):
        # displacement uniform in the ball of radius rmax
        while True:
            dv = (rng.random(3) - 0.5) * 2 * rmax
            if dv @ dv < rmax * rmax: break
        Ri, Rj = sim_rand_rot(2, rng)
        hs = []
        for q in prod_h:
            sh = shs[q]; Y = np.empty((sh['V'].shape[0], 3)); g = np.empty(sh['V'].shape[0])
            hs.append(bool(sim4.sheet_contact(Ri, dv, Rj, sh['V'], sh['E'], sh['stype'], sh['sp0'], sh['sp1'], sh['sp2'], sh['R0'],
                                              sh['xtab'], sh['utab'], Y, g, offs, 1.0, sr)))
        rotate_tris(T, Ri, np.zeros(3), TAr)
        rotate_tris(T * np.array([sr, sr, sr])[None, None, :] if sr != 1.0 else T, Rj, dv, TBr)
        hm = mesh_contact(TAr, TBr, cs, lo, nc)
        res.append((hs, bool(hm), float(np.linalg.norm(dv))))
    return res, rmax

if __name__ == "__main__":
    npairs = int(sys.argv[1]); seed = int(sys.argv[2]); out = []
    cases = [("flat", 0.0, 1.0, 1.0), ("a1.2/lam1.0", 1.2 / (2 * math.pi), 1.0, 1.0), ("P1,P2 a2.0/lam1.0", 2.0 / (2 * math.pi), 1.0, 1.0),
             ("P3 a1.5/lam0.7", 1.5 * 0.7 / (2 * math.pi), 0.7, 1.0), ("P4 a3.0/lam1.5", 3.0 * 1.5 / (2 * math.pi), 1.5, 1.0),
             ("a2.0/lam1.0, radius ratio 0.6", 2.0 / (2 * math.pi), 1.0, 0.6)]
    for name, amp, lam, sr in cases:
        for hm in (1 / 20, 1 / 40):
            hm = hm
            t0 = time.time(); res, rmax = run(amp, lam, hm, npairs, seed, sr)
            H = np.array([r[0] for r in res]); hmm = np.array([r[1] for r in res]); vol = 4 / 3 * math.pi * rmax ** 3
            row = dict(case=name, h_mesh=hm, npairs=npairs, contacts_mesh=int(hmm.sum()), Vex_mesh=float(hmm.mean() * vol),
                       Vex_mesh_err=float(vol * math.sqrt(hmm.mean() * (1 - hmm.mean()) / npairs)), secs=round(time.time() - t0))
            for c, q in enumerate((8, 16, 32)):
                row[f"prod_h{q}"] = dict(contacts=int(H[:, c].sum()), prod_only=int((H[:, c] & ~hmm).sum()), mesh_only=int((~H[:, c] & hmm).sum()),
                                         Vex=float(H[:, c].mean() * vol), ratio_to_mesh=float(H[:, c].mean() / hmm.mean()))
            print(json.dumps(row), flush=True); out.append(row)
    json.dump(out, open('ripple_mesh_test.json', 'w'), indent=1)
