"""Positive wrapping trials, translation invariance, face-contact convention, large-L disc means.
Run from review_work/: python test_followup.py
"""
import sys, os, json, math
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "code_copy", "literature_checks"))
from sim2e import disc_contact, random_rotations, run_once

rng = np.random.default_rng(7)

def supercell_Nc(pos, rot, rad, L):
    N = len(pos)
    images = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    bonds = []
    for i in range(N):
        ni = rot[i, :, 2]
        for j in range(i + 1, N):
            nj = rot[j, :, 2]
            for s in images:
                dv = pos[j] - pos[i] + np.array(s, float) * L
                if dv @ dv >= (rad[i] + rad[j]) ** 2:
                    continue
                if disc_contact(ni, float(rad[i]), dv, nj, float(rad[j])):
                    bonds.append((i, j, s))
    bonds.sort(key=lambda b: max(b[0], b[1]))
    parent = {}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def add(x):
        if x not in parent:
            parent[x] = x
    bi = 0
    for k in range(1, N + 1):
        add((k - 1, 0, 0, 0))
        while bi < len(bonds) and max(bonds[bi][0], bonds[bi][1]) < k:
            i, j, s = bonds[bi]; bi += 1
            for a, b, c in images:
                ja, jb, jc = a + s[0], b + s[1], c + s[2]
                if max(abs(ja), abs(jb), abs(jc)) > 1:
                    continue
                u = (i, a, b, c); v = (j, ja, jb, jc)
                add(u); add(v)
                ru, rv = find(u), find(v)
                if ru != rv:
                    parent[ru] = rv
        for i in range(k):
            home = find((i, 0, 0, 0))
            for s in images:
                if s == (0, 0, 0):
                    continue
                node = (i, s[0], s[1], s[2])
                if node in parent and find(node) == home:
                    return k
    return -1

print("positive wraps, L=7, N=120")
mism = 0
base_nc = None
pos0 = rot0 = None
for trial in range(4):
    L = 7.0
    n = 120
    pos = rng.random((n, 3)) * L
    rot = random_rotations(n, rng)
    rad = np.ones(n)
    a = int(run_once(0, pos, rot, rad, np.zeros(n), L, rad))
    b = supercell_Nc(pos, rot, rad, L)
    print(f"  trial {trial}: sim {a} supercell {b} {'OK' if a == b else 'DIFF'}")
    mism += a != b
    if base_nc is None and a > 0:
        base_nc, pos0, rot0 = a, pos, rot
print("mismatches", mism)

if pos0 is not None:
    L = 7.0
    rad = np.ones(len(pos0))
    shift = (pos0 + np.array([0.37, 1.2, 2.5])) % L
    nc_s = int(run_once(0, shift, rot0, rad, np.zeros(len(pos0)), L, rad))
    box = (pos0 + np.array([L, 0, 0]))
    nc_b = int(run_once(0, box, rot0, rad, np.zeros(len(pos0)), L, rad))
    Rz = np.array([[0., -1, 0], [1., 0, 0], [0, 0, 1.]])
    pos_r = (pos0 @ Rz.T) % L
    rot_r = np.einsum("ij,njk->nik", Rz, rot0)
    nc_r = int(run_once(0, pos_r, rot_r, rad, np.zeros(len(pos0)), L, rad))
    print(f"invariance on a percolating draw: base {base_nc} shift {nc_s} +L {nc_b} rot90 {nc_r}")

# face extent versus chord-in-square
print("\nface convention")
def chord_hits_square(c, n, R, face_coord, axis, L):
    """Does the disc intersect the square face sitting at coordinate face_coord on `axis`?"""
    e = np.zeros(3); e[axis] = 1.0
    # signed distance from centre to the plane
    dist = (c[axis] - face_coord)
    # intersection with the infinite plane
    # direction u = n × e
    u = np.cross(n, e)
    us = np.linalg.norm(u)
    if us < 1e-12:
        # disc parallel to the face
        if abs(dist) > 1e-10:
            return False
        # disc lies in the plane: circle-square overlap
        axes = [q for q in range(3) if q != axis]
        # sample the circle; exact SAT in 2D
        cc = c[axes]
        # square [0,L]^2
        # circle overlaps square if closest point of square to centre is within R
        closest = np.clip(cc, 0, L)
        return np.sum((cc - closest) ** 2) <= R * R + 1e-10
    u = u / us
    # point on the line: solve (x-c)·n = 0 and x[axis]=face_coord, x·u free
    # particular solution: move from c along the plane of the disc toward the face plane
    # x = c + a * (e - (e·n) n)  chosen so x[axis]=face_coord, and stay in disc plane
    en = e - (e @ n) * n
    # en[axis] = 1 - n[axis]^2
    denom = en[axis]
    if abs(denom) < 1e-14:
        return False
    a = (face_coord - c[axis]) / denom
    p = c + a * en  # on both planes, and p is the point of the line closest to c? not always, but it is on the line
    # half chord
    # distance from c to the line. p should be the projection if en is in the disc plane and we solved the plane.
    # Verify (p-c)·n ~ 0
    d2 = np.sum((p - c) ** 2)
    if d2 > R * R + 1e-8:
        return False
    half = math.sqrt(max(0.0, R * R - d2))
    # segment p ± half u, in the two in-plane coordinates, versus the square
    axes = [q for q in range(3) if q != axis]
    def inside_seg():
        # overlap of segment with the square: sample densely plus endpoints, and also clip
        ts = np.linspace(-half, half, 41)
        for t in ts:
            q = p + t * u
            if 0 <= q[axes[0]] <= L and 0 <= q[axes[1]] <= L:
                return True
        return False
    return inside_seg()

def extent_hits(c, n, R, L):
    hits = []
    for a in range(3):
        e = R * math.sqrt(max(0.0, 1 - n[a] ** 2))
        lo = c[a] - e <= 0
        hi = c[a] + e >= L
        hits.append((lo, hi))
    return hits

L = 16.0
N = 4000
pos = rng.random((N, 3)) * L
rot = random_rotations(N, rng)
extent_only = 0; chord_only = 0; both = 0; extent_n = 0
for i in range(N):
    nrm = rot[i, :, 2]; c = pos[i]
    ext = extent_hits(c, nrm, 1.0, L)
    for a in range(3):
        for side, coord in ((0, 0.0), (1, L)):
            if not ext[a][side]:
                continue
            extent_n += 1
            ok = chord_hits_square(c, nrm, 1.0, coord, a, L)
            if ok:
                both += 1
            else:
                extent_only += 1
print(f"extent face marks {extent_n}, of which chord misses {extent_only} ({100*extent_only/max(extent_n,1):.2f}%)")

# large-box disc means in doc 02
print("\nlarge L disc files")
base = os.path.join(ROOT, "code_copy", "verification")
acc = []
thin = json.load(open(os.path.join(base, "fss_thin.json")))
for k, v in thin.items():
    a = np.array(v, float)
    m = a.mean(); se = a.std(ddof=1) / math.sqrt(len(a))
    print(f"  fss_thin.json L={k}: {m:.4f} +- {se:.4f} n={len(a)}")
    acc.append((float(k), m, se))
for L, fn in ((64, "fss_L64.json"), (96, "fss_L96.json"), (128, "fss_L128.json")):
    a = np.array(json.load(open(os.path.join(base, fn))), float)
    m = a.mean(); se = a.std(ddof=1) / math.sqrt(len(a))
    print(f"  {fn} L={L}: {m:.4f} +- {se:.4f} n={len(a)}")
    acc.append((float(L), m, se))
for label, pred in (("L>=32", lambda L: L >= 32), ("four largest", None), ("L>=48", lambda L: L >= 48)):
    rows = sorted(acc, key=lambda t: t[0])
    if label == "four largest":
        rows = rows[-4:]
    else:
        rows = [t for t in rows if pred(t[0])]
    w = np.array([1 / se ** 2 for _, _, se in rows]); y = np.array([m for _, m, _ in rows])
    print(f"  weighted {label}: {np.sum(w*y)/np.sum(w):.4f} +- {1/math.sqrt(w.sum()):.4f} from {[t[0] for t in rows]}")

# sphere L=32
sph = json.load(open(os.path.join(base, "fss_spheres.json")))
for L in ("16.0", "24.0", "32.0"):
    a = np.array(sph[L]); print(f"sphere exact L={L}: {a.mean():.4f} +- {a.std(ddof=1)/math.sqrt(len(a)):.4f} n={len(a)}")

# F6 bounding radii versus L/2
print("\nF6/F5 bounding radius versus L/2")
sys.path.insert(0, os.path.join(ROOT, "code_copy"))
os.chdir(os.path.join(ROOT, "code_copy"))
import exam_sim
for case, sig in (("F5", None), ("F6", None)):
    rng2 = np.random.default_rng(1)
    R, t, kap, ell = exam_sim.build(case, 200000, rng2)
    hth = t / 2
    brad = np.sqrt(R ** 2 + hth ** 2)
    print(f"  {case}: n={len(R)} maxR={R.max():.3f} max t={t.max():.3f} max brad={brad.max():.3f} "
          f"2*maxbrad={2*brad.max():.3f} L/2=16  margin L/2-2*max={16-2*brad.max():.3f}")
