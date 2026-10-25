"""Adversarial geometry tests: cell list, winding, GJK, discs, polygons, spheroids, face contact.
Run from review_work/: python test_geometry.py
"""
import sys, os, math
import numpy as np
from scipy.optimize import minimize

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "code_copy", "literature_checks"))
sys.path.insert(0, os.path.join(ROOT, "code_copy"))

from sim2e import (disc_contact, gjk_cyl, gjk_ell, sup_cyl, sup_ell, poly_contact,
                   random_rotations, first_perc, run_once, EMPTY)
from poly_test import polygon

def section(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72, flush=True)

rng = np.random.default_rng(20261025)

# ---------------------------------------------------------------- cell list
section("cell list versus all pairs")

def cell_pairs(pos, brad, L, ncell):
    """Pairs (j < i) the neighbour loop of first_perc would test, using the same cell rule."""
    cs = L / ncell
    N = len(pos)
    head = {(a, b, c): [] for a in range(ncell) for b in range(ncell) for c in range(ncell)}
    seen = []
    for i in range(N):
        cx, cy, cz = [int(pos[i, q] / cs) % ncell for q in range(3)]
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in head[((cx + dx) % ncell, (cy + dy) % ncell, (cz + dz) % ncell)]:
                        dv0 = pos[j] - pos[i]
                        sh = -np.round(dv0 / L)
                        dv = dv0 + sh * L
                        if dv @ dv < (brad[i] + brad[j]) ** 2:
                            seen.append((j, i))
        head[(cx, cy, cz)].append(i)
    return set(seen)

def all_pairs(pos, brad, L):
    N = len(pos)
    out = set()
    for i in range(N):
        dv0 = pos[:i] - pos[i]
        sh = -np.round(dv0 / L)
        dv = dv0 + sh * L
        rr = brad[:i] + brad[i]
        hit = np.where((dv * dv).sum(1) < rr * rr)[0]
        for j in hit:
            out.add((int(j), i))
    return out

def check_cells(label, pos, brad, L):
    bmax = float(brad.max())
    ncell = max(3, int(L // (2 * bmax)))
    cs = L / ncell
    A = cell_pairs(pos, brad, L, ncell)
    B = all_pairs(pos, brad, L)
    miss = B - A
    extra = A - B
    # T5-style vacuity: particles with brad > cs - bmax
    n_big = int(np.sum(brad > (cs - bmax)))
    print(f"{label}: N={len(pos)} ncell={ncell} cs={cs:.3f} 2*bmax={2*bmax:.3f} "
          f"pairs {len(B)} missed {len(miss)} extra {len(extra)}  T5-big-particles {n_big}")
    return len(miss)

# monodisperse L=32, but only 400 particles so all-pairs is cheap; cell rule still ncell>3
pos = rng.random((400, 3)) * 32
brad = np.full(400, math.sqrt(1 + 0.02 ** 2))
check_cells("mono L32 N400", pos, brad, 32.0)

# broad radii, ncell>3: L=48, R up to 3, 2*bmax~6.1, floor(48/6.1)=7
R = rng.lognormal(-0.5 * 0.45 ** 2, 0.45, size=500)
R = np.clip(R, 0.15, 3.0)
brad = np.sqrt(R ** 2 + 0.05 ** 2)
pos = rng.random((500, 3)) * 48
check_cells("broad L48 cap3", pos, brad, 48.0)

# particles jammed on a cell boundary
L = 30.0
bmax = 1.2
ncell = max(3, int(L // (2 * bmax)))
cs = L / ncell
pos = rng.random((300, 3)) * L
# put 30 particles exactly on planes k*cs
for k in range(30):
    pos[k, 0] = (k % ncell) * cs
brad = rng.uniform(0.3, bmax, size=300)
check_cells("boundary L30", pos, brad, L)

# minimum-image double contact (outside the safe regime)
section("minimum-image double contact, two large discs")
L = 10.0
# two coplanar discs, R=4, centres 3 apart: both the near image and the far image overlap
pos = np.array([[0.0, 5.0, 5.0], [3.0, 5.0, 5.0]])
rad = np.array([4.0, 4.0])
hth = np.zeros(2)
rot = np.zeros((2, 3, 3)); rot[:, 2, 2] = 1.0; rot[:, 0, 0] = 1.0; rot[:, 1, 1] = 1.0
brad = rad.copy()
nc = run_once(0, pos, rot, rad, hth, L, brad)
# explicit second image
d_near = 3.0
d_far = L - 3.0
print(f"near {d_near} far {d_far} Rsum {8.0}  first_perc -> {nc} ( -1 means no wrap detected )")
print(f"disc near contact {disc_contact(np.array([0.,0.,1.]), 4., np.array([3.,0.,0.]), np.array([0.,0.,1.]), 4.)}")
print(f"disc far  contact {disc_contact(np.array([0.,0.,1.]), 4., np.array([3.-L,0.,0.]), np.array([0.,0.,1.]), 4.)}")

# ---------------------------------------------------------------- winding versus 3x3x3 supercell
section("winding versus supercell")

def supercell_Nc(pos, rot, rad, L):
    """First k such that the contact graph on a 3x3x3 supercell joins a particle to its own nonzero image.
    Disc mode. Independent of the offset union-find: BFS on explicit copies."""
    N = len(pos)
    images = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    # precompute, for every pair i<j, the image shifts of j relative to i that produce a contact
    # shift s means pos[j]+s*L contacts pos[i]
    bonds = []  # (i, j, sx, sy, sz) with i<j, contact between i at 0 and j at s
    for i in range(N):
        ni = rot[i, :, 2]
        for j in range(i + 1, N):
            nj = rot[j, :, 2]
            for s in images:
                dv = pos[j] - pos[i] + np.array(s) * L
                if dv @ dv >= (rad[i] + rad[j] + 1e-9) ** 2:
                    continue
                if disc_contact(ni, rad[i], dv, nj, rad[j]):
                    bonds.append((i, j, s[0], s[1], s[2]))
    # incremental: a bond is active when both endpoints have been inserted (index < k, since 0..k-1)
    # We add bonds in order of max(i,j)
    bonds.sort(key=lambda b: max(b[0], b[1]))
    parent = {}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def add_node(x):
        if x not in parent:
            parent[x] = x
    bi = 0
    for k in range(1, N + 1):
        add_node((k - 1, 0, 0, 0))
        while bi < len(bonds) and max(bonds[bi][0], bonds[bi][1]) < k:
            i, j, sx, sy, sz = bonds[bi]; bi += 1
            # bond connects every copy (i, a) to (j, a+s) for images we keep inside -1..1
            for a, b, c in images:
                ja, jb, jc = a + sx, b + sy, c + sz
                if max(abs(ja), abs(jb), abs(jc)) > 1:
                    continue
                u = (i, a, b, c); v = (j, ja, jb, jc)
                add_node(u); add_node(v)
                ru, rv = find(u), find(v)
                if ru != rv:
                    parent[ru] = rv
        # wrap if some inserted particle's home copy meets a nonzero image
        for i in range(k):
            home = find((i, 0, 0, 0))
            for s in images:
                if s == (0, 0, 0):
                    continue
                node = (i, s[0], s[1], s[2])
                if node in parent and find(node) == home:
                    return k
    return -1

mism = 0
for trial in range(8):
    L = 8.0
    n = 70
    pos = rng.random((n, 3)) * L
    rot = random_rotations(n, rng)
    rad = np.ones(n)
    hth = np.zeros(n)
    a = int(run_once(0, pos, rot, rad, hth, L, rad))
    b = supercell_Nc(pos, rot, rad, L)
    flag = "OK" if a == b else "DIFF"
    if a != b:
        mism += 1
    print(f"  trial {trial}: sim {a} supercell {b} {flag}")
print(f"mismatches {mism}/8")

# crafted: chain that closes with two different offsets
section("crafted double-offset wrap")
L = 20.0
# A at x=2, B at x=6 (shift 0), C at x=17 contacts A across the boundary, D at x=11.5 contacts B and C
pts = np.array([[2.0, 10, 10], [6.0, 10, 10], [17.0, 10, 10], [11.5, 10, 10]], float)
# radii 3 so distances 4, 5.5, 5.5 work and the boundary gap from 17 to 2+20=22 is 5
rad = np.full(4, 3.0)
rot = np.zeros((4, 3, 3)); rot[:, 2, 2] = 1; rot[:, 0, 0] = 1; rot[:, 1, 1] = 1
hth = np.zeros(4)
a = int(run_once(0, pts, rot, rad, hth, L, rad))
b = supercell_Nc(pts, rot, rad, L)
print(f"crafted chain sim {a} supercell {b}")
# also the two-particle double-image case through the supercell (R=4, L=10)
pos2 = np.array([[1.0, 5, 5], [4.0, 5, 5]], float)
rad2 = np.full(2, 4.0)
rot2 = np.zeros((2, 3, 3)); rot2[:, 2, 2] = 1; rot2[:, 0, 0] = 1; rot2[:, 1, 1] = 1
print(f"two large discs sim {int(run_once(0, pos2, rot2, rad2, np.zeros(2), 10.0, rad2))} "
      f"supercell {supercell_Nc(pos2, rot2, rad2, 10.0)}")

# translation / 90-degree rotation invariance on one random config
section("invariance")
L = 9.0
n = 80
pos = rng.random((n, 3)) * L
rot = random_rotations(n, rng)
rad = np.ones(n); hth = np.zeros(n)
base = int(run_once(0, pos, rot, rad, hth, L, rad))
pos_t = (pos + np.array([L, 0, 0])) % (2 * L)  # shift by one box, leave in [0, 2L)
# run_once cell index uses int(pos/cs) % ncell, so [L, 2L) should map like [0, L)
nc_t = int(run_once(0, (pos + np.array([0.3, 0.7, 1.1])) % L, rot, rad, hth, L, rad))
# 90-degree rotation about z: (x,y)->(-y,x), and rotate the normal the same way
Rz = np.array([[0., -1, 0], [1., 0, 0], [0, 0, 1]])
pos_r = pos @ Rz.T
rot_r = np.einsum("ij,njk->nik", Rz, rot)
# keep inside box
pos_r = pos_r % L
nc_r = int(run_once(0, pos_r, rot_r, rad, hth, L, rad))
# swap argument symmetry of disc_contact on 2000 pairs
bad_swap = 0
for _ in range(2000):
    d = rng.normal(size=3); n1 = rng.normal(size=3); n2 = rng.normal(size=3)
    n1 /= np.linalg.norm(n1); n2 /= np.linalg.norm(n2)
    r1, r2 = rng.uniform(0.4, 1.6, size=2)
    A = disc_contact(n1, r1, d, n2, r2)
    B = disc_contact(n2, r2, -d, n1, r1)
    bad_swap += (A != B)
print(f"base Nc {base}  translate {nc_t}  rot90 {nc_r}  disc swap mismatches {bad_swap}/2000")

# order dependence (expected): reverse the insertion order
pos_rev = pos[::-1].copy(); rot_rev = rot[::-1].copy()
nc_rev = int(run_once(0, pos_rev, rot_rev, rad, hth, L, rad))
print(f"reversed insertion order Nc {nc_rev} (base {base}); difference is the estimator, not a contact change")

# ---------------------------------------------------------------- disc degenerates
section("disc degenerates")
ez = np.array([0., 0., 1.])
ex = np.array([1., 0., 0.])
# coplanar overlap, edge touch, separated, parallel separated, nearly parallel intersecting
cases = []
cases.append(("coplanar overlap", disc_contact(ez, 1, np.array([0.5, 0, 0.]), ez, 1), True))
cases.append(("coplanar edge", disc_contact(ez, 1, np.array([2.0, 0, 0.]), ez, 1), False))  # strict <
cases.append(("coplanar just inside", disc_contact(ez, 1, np.array([1.999, 0, 0.]), ez, 1), True))
cases.append(("parallel gap", disc_contact(ez, 1, np.array([0.0, 0, 0.1]), ez, 1), False))
cases.append(("antiparallel overlap", disc_contact(ez, 1, np.array([0.5, 0, 0.]), -ez, 1), True))
# perpendicular discs, centres separated by 0.5 along x: planes x=0? disc i normal z through origin, disc j normal x at (0.5,0,0)
# intersection line is parallel to y through (0.5, 0, 0) wait plane j: (X - d)·ex = 0 => X_x = 0.5
# plane i: Z=0. Line x=0.5, z=0. Distance from i centre 0.5 < 1, from j centre 0 < 1. Overlap.
cases.append(("perpendicular overlap", disc_contact(ez, 1, np.array([0.5, 0, 0.]), ex, 1), True))
# nearly parallel, actually intersecting: angle 1e-6, separation of planes small enough
n2 = np.array([math.sin(1e-6), 0, math.cos(1e-6)])
cases.append(("near-parallel intersect", disc_contact(ez, 1, np.array([0.2, 0, 1e-8]), n2, 1), True))
# angle below the 1e-9 cutoff, coplanar overlap should still be caught by the parallel branch
n3 = np.array([math.sin(1e-12), 0, math.cos(1e-12)])
cases.append(("below cutoff coplanar overlap", disc_contact(ez, 1, np.array([0.2, 0, 0.0]), n3, 1), True))
# below cutoff but planes separated by 1e-8: parallel branch rejects (h>1e-12) even if a tiny angle would still intersect
cases.append(("below cutoff tiny gap", disc_contact(ez, 1, np.array([0.2, 0, 1e-8]), n3, 1), False))
for name, got, expect in cases:
    print(f"  {name}: got {got} expect {expect} {'OK' if got == expect else 'DIFF'}")

# ---------------------------------------------------------------- GJK spheres, cylinders, symmetry, support
section("GJK sphere (R=h) versus centre distance")
bad = 0; n = 3000; gaps = []
for i in range(n):
    c = rng.normal(size=3) * rng.uniform(0.2, 3.0)
    # random radii
    R1 = rng.uniform(0.3, 1.5); R2 = rng.uniform(0.3, 1.5)
    n1 = rng.normal(size=3); n1 /= np.linalg.norm(n1)
    n2 = rng.normal(size=3); n2 /= np.linalg.norm(n2)
    dist = np.linalg.norm(c)
    truth = dist < R1 + R2
    pred = gjk_ell(np.zeros(3), n1, R1, R1, c, n2, R2, R2)
    if pred != truth:
        bad += 1
        gaps.append(dist - (R1 + R2))
print(f"sphere mismatches {bad}/{n}" + (f"  gap range of mismatches {min(gaps):.3e} .. {max(gaps):.3e}" if gaps else ""))

section("GJK cylinder support and overlap oracle")

def support_sample(kind, n, R, h, direction, ntheta=400, nax=80):
    d = np.asarray(direction, float)
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    # orthonormal frame
    a = np.array([1., 0, 0]) if abs(n[0]) < 0.9 else np.array([0., 1., 0])
    e1 = np.cross(n, a); e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1)
    best = -1e300
    th = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    if kind == "cyl":
        for z in np.linspace(-h, h, nax):
            for t in th:
                p = z * n + R * (math.cos(t) * e1 + math.sin(t) * e2)
                best = max(best, float(p @ d))
    else:
        for zf in np.linspace(-1, 1, nax):
            rad = R * math.sqrt(max(0.0, 1 - zf * zf))
            z = h * zf
            for t in th:
                p = z * n + rad * (math.cos(t) * e1 + math.sin(t) * e2)
                best = max(best, float(p @ d))
    return best

def sup_value_cyl(n, R, h, d):
    out = np.empty(3)
    sup_cyl(np.zeros(3), np.asarray(n, float), R, h, np.asarray(d, float), out)
    return float(out @ np.asarray(d, float))

sup_bad = 0
dirs = [np.array([0., 0., 1.]), np.array([1., 0., 0.]), np.array([1., 1., 1.]),
        np.array([1e-8, 0, 1.]), np.array([1., 0, 1e-8])]
for h in (0.05, 0.5, 2.0):
    for d in dirs:
        nrm = np.array([0., 0., 1.])
        analytic = sup_value_cyl(nrm, 1.0, h, d)
        samp = support_sample("cyl", nrm, 1.0, h, d)
        rel = abs(analytic - samp) / (abs(samp) + 1e-12)
        ok = rel < 1e-3 or abs(analytic - samp) < 1e-4
        sup_bad += not ok
        if not ok:
            print(f"  SUPPORT DIFF h={h} d={d} analytic {analytic:.6f} sample {samp:.6f}")
print(f"cylinder support disagreements {sup_bad}")

def cyl_gap(cB, nA, RA, hA, nB, RB, hB):
    """Minimum distance between two solid cylinders; negative means overlap depth (approx 0 if overlapping, we return signed).
    Convex SLSQP. Returns distance (0 if overlapping within tolerance)."""
    nA = np.asarray(nA, float); nB = np.asarray(nB, float)
    cB = np.asarray(cB, float)
    def obj(v):
        return np.sum((v[:3] - v[3:]) ** 2)
    cons = []
    def pack(c, n, R, h, sl):
        def axial(v):
            x = v[sl:sl + 3]
            return h - abs(n @ (x - c))
        def radial(v):
            x = v[sl:sl + 3]
            q = (x - c) - n * (n @ (x - c))
            return R ** 2 - q @ q
        cons.append({"type": "ineq", "fun": axial})
        cons.append({"type": "ineq", "fun": radial})
    pack(np.zeros(3), nA, RA, hA, 0)
    pack(cB, nB, RB, hB, 3)
    x0 = np.concatenate([np.zeros(3), cB])
    res = minimize(obj, x0, method="SLSQP", constraints=cons, tol=1e-12,
                   options={"maxiter": 200, "ftol": 1e-14})
    if not res.success:
        return None
    return math.sqrt(max(res.fun, 0.0))

# controlled side-side gaps, parallel axes
print("parallel cylinders, gap from -0.05 to 0.05 (negative = overlap)")
for h in (1e-3, 0.2, 2.0):
    for gap in (-0.05, -0.01, -1e-4, 0.0, 1e-4, 0.01, 0.05):
        cB = np.array([2 + gap, 0.0, 0.0])  # radii 1+1
        g = gjk_cyl(np.zeros(3), np.array([0., 0., 1.]), 1.0, h, cB, np.array([0., 0., 1.]), 1.0, h)
        print(f"  h={h:<6} gap={gap:+.4f} gjk={g}")

# random pairs versus SLSQP, skip the tangent band
bad_g = 0; unsure = 0; tested = 0
swap_bad = 0
for i in range(200):
    n1 = rng.normal(size=3); n1 /= np.linalg.norm(n1)
    n2 = rng.normal(size=3); n2 /= np.linalg.norm(n2)
    c = rng.normal(size=3); c *= rng.uniform(0.3, 2.5) / (np.linalg.norm(c) + 1e-15)
    h1 = rng.choice([1e-3, 0.05, 0.5, 2.0]); h2 = rng.choice([1e-3, 0.05, 0.5, 2.0])
    R1 = rng.uniform(0.5, 1.2); R2 = rng.uniform(0.5, 1.2)
    dist = cyl_gap(c, n1, R1, h1, n2, R2, h2)
    pred = bool(gjk_cyl(np.zeros(3), n1, R1, h1, c, n2, R2, h2))
    pred2 = bool(gjk_cyl(c, n2, R2, h2, np.zeros(3), n1, R1, h1))
    swap_bad += pred != pred2
    if dist is None:
        unsure += 1
        continue
    tested += 1
    if dist > 2e-3 and pred:
        bad_g += 1
        print(f"  FALSE contact dist {dist:.4f} h {h1},{h2}")
    if dist < -1:  # unreachable, distance is >= 0
        pass
    # overlap: if the minimizer found a positive distance, they are separate.
    # If distance is ~0 they may be overlapping OR tangent. Check a shrunk pair.
    if dist <= 1e-4:
        # bodies reach each other. Shrink radii by 0.01 and see if still contacting analytically via another minimization
        dist2 = cyl_gap(c, n1, max(R1 - 0.02, 0.05), h1, n2, max(R2 - 0.02, 0.05), h2)
        if dist2 is not None and dist2 > 5e-3 and pred:
            # original was a near touch that GJK called contact; acceptable
            pass
        if dist2 is not None and dist2 <= 1e-4 and not pred:
            bad_g += 1
            print(f"  MISSED overlap (still touching after 0.02 shrink) h {h1} {h2}")
print(f"cylinder oracle mismatches {bad_g} in {tested} (optimizer failed {unsure}); swap mismatches {swap_bad}/200")

# ---------------------------------------------------------------- spheroid h=0.001
section("spheroid h=0.001 and support")
def ell_support_sample(n, R, h, d):
    return support_sample("ell", n, R, h, d)

ell_sup_bad = 0
for h in (1.0, 0.1, 0.001):
    for d in (np.array([0., 0., 1.]), np.array([1., 0., 0.]), np.array([0.2, 0.1, 1.]), np.array([1., 0., 1e-6])):
        nrm = np.array([0., 0., 1.])
        out = np.empty(3)
        sup_ell(np.zeros(3), nrm, 1.0, h, d.astype(float), out)
        analytic = float(out @ d)
        samp = ell_support_sample(nrm, 1.0, h, d)
        if abs(analytic - samp) > 1e-3 * (abs(samp) + 1):
            ell_sup_bad += 1
            print(f"  ELL SUPPORT h={h} d={d} analytic {analytic:.6f} sample {samp:.6f}")
print(f"spheroid support disagreements {ell_sup_bad}")

# versus disc contact
bad_d = 0; n = 1500
for i in range(n):
    d = rng.normal(size=3); d *= rng.uniform(0.0, 2.2) / (np.linalg.norm(d) + 1e-15)
    n1 = rng.normal(size=3); n1 /= np.linalg.norm(n1)
    n2 = rng.normal(size=3); n2 /= np.linalg.norm(n2)
    A = bool(disc_contact(n1, 1.0, d, n2, 1.0))
    B = bool(gjk_ell(np.zeros(3), n1, 1.0, 0.001, d, n2, 1.0, 0.001))
    if A != B:
        bad_d += 1
print(f"spheroid h=0.001 vs disc mismatches {bad_d}/{n}")

# a known separated face-on pair and an overlapping edge-on pair
nA = np.array([0., 0., 1.])
print("face-on centres 0.01 apart along normal, h=0.001 (overlap expected)",
      bool(gjk_ell(np.zeros(3), nA, 1.0, 0.001, np.array([0, 0, 0.0015]), nA, 1.0, 0.001)))
print("face-on centres 0.01 apart (separated)",
      bool(gjk_ell(np.zeros(3), nA, 1.0, 0.001, np.array([0, 0, 0.01]), nA, 1.0, 0.001)))
print("edge-on centres 1.5 apart (overlap)",
      bool(gjk_ell(np.zeros(3), nA, 1.0, 0.001, np.array([1.5, 0, 0]), nA, 1.0, 0.001)))
print("edge-on centres 2.1 apart (separated)",
      bool(gjk_ell(np.zeros(3), nA, 1.0, 0.001, np.array([2.1, 0, 0]), nA, 1.0, 0.001)))

# ---------------------------------------------------------------- polygons
section("polygons")
Vsq = polygon("square"); Vtr = polygon("triangle")
print("square area", 0.5 * abs(np.sum(Vsq[:, 0] * np.roll(Vsq[:, 1], -1) - np.roll(Vsq[:, 0], -1) * Vsq[:, 1])))
print("triangle area", 0.5 * abs(np.sum(Vtr[:, 0] * np.roll(Vtr[:, 1], -1) - np.roll(Vtr[:, 0], -1) * Vtr[:, 1])))

def world(Rm, c, V):
    return V @ Rm.T + c  # Rm columns are local axes, as in the code (world = Rm @ local)

def my_interval(Rm, c, V, pn, pc, u):
    X = world(Rm, c, V)
    s = (X - pc) @ pn
    ts = []
    nv = len(V)
    for k in range(nv):
        k2 = (k + 1) % nv
        # count a crossing when the open segment crosses or touches from either side, once per vertex on the plane
        if s[k] == 0 and s[k2] == 0:
            ts.append(float(X[k] @ u)); ts.append(float(X[k2] @ u))
            continue
        if (s[k] < 0 and s[k2] > 0) or (s[k] > 0 and s[k2] < 0) or s[k] == 0 or s[k2] == 0:
            if s[k] == 0 and s[k2] == 0:
                continue
            # if only the endpoint is on the plane, record it once (when k is the on-plane vertex)
            if s[k] == 0:
                ts.append(float(X[k] @ u)); continue
            if s[k2] == 0:
                continue  # recorded when that vertex is k
            f = s[k] / (s[k] - s[k2])
            p = X[k] + f * (X[k2] - X[k])
            ts.append(float(p @ u))
    if len(ts) < 2:
        return None
    return min(ts), max(ts)

def my_contact(Ri, dv, Rj, V):
    ni = Ri[:, 2]; nj = Rj[:, 2]
    u = np.cross(ni, nj); us = np.linalg.norm(u)
    if us < 1e-12:
        # coplanar test: 2D separating axis on the common plane, only if plane separation is ~0
        if abs(np.dot(dv, ni)) > 1e-9:
            return False
        Xi = world(Ri, np.zeros(3), V)[:, :2] if False else None
        # project both polygons onto the plane basis
        a = np.array([1., 0, 0]) if abs(ni[0]) < 0.9 else np.array([0., 1., 0])
        e1 = np.cross(ni, a); e1 /= np.linalg.norm(e1); e2 = np.cross(ni, e1)
        def proj(Rm, c):
            X = world(Rm, c, V)
            return np.column_stack([X @ e1, X @ e2])
        P = proj(Ri, np.zeros(3)); Q = proj(Rj, dv)
        # SAT
        def axes(poly):
            d = np.roll(poly, -1, axis=0) - poly
            return np.column_stack([-d[:, 1], d[:, 0]])
        for ax in np.vstack([axes(P), axes(Q)]):
            na = np.linalg.norm(ax)
            if na < 1e-15:
                continue
            ax = ax / na
            p0, p1 = (P @ ax).min(), (P @ ax).max()
            q0, q1 = (Q @ ax).min(), (Q @ ax).max()
            if p1 < q0 - 1e-10 or q1 < p0 - 1e-10:
                return False
        return True
    u = u / us
    a = my_interval(Ri, np.zeros(3), V, nj, dv, u)
    b = my_interval(Rj, dv, V, ni, np.zeros(3), u)
    if a is None or b is None:
        return False
    return max(a[0], b[0]) <= min(a[1], b[1]) + 1e-9

# random agreement
bad_p = 0; n = 2000; both = 0
rots = random_rotations(2 * n, rng)
dvs = rng.normal(size=(n, 3))
for i in range(n):
    A = bool(poly_contact(rots[2 * i], dvs[i], rots[2 * i + 1], Vsq))
    B = bool(my_contact(rots[2 * i], dvs[i], rots[2 * i + 1], Vsq))
    both += A
    if A != B:
        bad_p += 1
print(f"square random vs independent: mismatches {bad_p}/{n} ({both} code contacts)")

# crafted degenerates
I = np.eye(3)
# coplanar overlapping squares
print("coplanar overlap code", bool(poly_contact(I, np.array([0.2, 0.1, 0.0]), I, Vsq)),
      "independent", my_contact(I, np.array([0.2, 0.1, 0.0]), I, Vsq))
# vertex of j in the plane of i, rest of j strictly on +z (positive side of ni=+z)
# Build Rj so local vertex 0 goes to a point in z=0 and the polygon tilts upward.
# Simpler: put j's vertices by choosing a rotation whose local xy-plane is tilted.
# Local vertex v0 of the square. Use a rotation by a small angle about y: Ry
# Generic rotation so a single vertex is the unique extreme.
ax = np.array([0.3, 0.5, 0.8]); ax /= np.linalg.norm(ax)
ang = 0.9; K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
Ry = np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * (K @ K)
Vj = world(Ry, np.zeros(3), Vsq)
# Put the unique extreme vertex exactly at the origin (inside the other square).
kmin = int(np.argmin(Vj[:, 2])); dv = -Vj[kmin]
Zp = Vj + dv
print("vertex-in-plane, polygon on + side: code", bool(poly_contact(I, dv, Ry, Vsq)),
      "indep", my_contact(I, dv, Ry, Vsq))
print("   positions", np.round(Zp, 4))
kmax = int(np.argmax(Vj[:, 2])); dvn = -Vj[kmax]
Zn = Vj + dvn
print("vertex-in-plane, polygon on - side: code", bool(poly_contact(I, dvn, Ry, Vsq)),
      "indep", my_contact(I, dvn, Ry, Vsq))
print("   positions", np.round(Zn, 4))

# edge in the plane: rotate 90 deg about y so the square stands in the xz plane? 
# Take two vertices at the same z by a  rotation about x of the square, then shift.
Rx = np.array([[1, 0, 0], [0, math.cos(0.7), -math.sin(0.7)], [0, math.sin(0.7), math.cos(0.7)]])
Ve = world(Rx, np.zeros(3), Vsq)
# an edge of the original square has two vertices; after Rx their z differ.
# Force an edge into the plane by rotating about the edge. Use a house-built polygon: 
# identity square, and a second square rotated 90 deg about x (stands vertical) shifted so its bottom edge is in z=0.
Rv = np.array([[1., 0, 0], [0, 0, -1], [0, 1, 0]])  # y->z, z->-y? columns: local x->x, local y->z, local z->-y
# world = Rv @ local. local y becomes world z. Square vertices (±h,±h,0): z_world = local y component through column 1.
# Rv[:,1] = (0,0,1), so world z = local y. Vertices with local y = -h have world z = -h.
h = abs(Vsq[0, 0])
dv_e = np.array([0.15, 0.0, h])  # lift by h so the low edge (local y=-h) sits at z=0
Ze = world(Rv, dv_e, Vsq)
print("edge-in-plane z values", np.round(Ze[:, 2], 5))
print("edge-in-plane code", bool(poly_contact(I, dv_e, Rv, Vsq)), "indep", my_contact(I, dv_e, Rv, Vsq))
# flip normal of the vertical square so the body is on the other side of the plane... 
# The vertical square has points at z>=0 after the shift. To put the body on z<=0, shift the opposite way
# and flip: use -shift of the high edge.
dv_em = np.array([0.15, 0.0, -h])
Ze2 = world(Rv, dv_em, Vsq)
print("edge other side z", np.round(Ze2[:, 2], 5), "code", bool(poly_contact(I, dv_em, Rv, Vsq)),
      "indep", my_contact(I, dv_em, Rv, Vsq))

print("\nDONE")
