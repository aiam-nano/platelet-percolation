"""Degenerate and random cases of the zero-thickness polygon contact (sim2e.poly_contact, mode 4).

Added in v0.1.1 after the Level 4 review (finding F3).  Checks, against an independent oracle, that
  (1) exactly parallel, coplanar plates in contact are detected, separated or merely touching ones are not;
  (2) a vertex or an edge lying in the other plate's plane gives the same answer from either side and with the
      roles of the two plates swapped;
  (3) random non-parallel pairs agree with an independent chord-overlap test.
Run from the repository root:  python tests/test_polygon_contact.py     (about 20 s, exit status 1 on failure)
"""
import os, sys, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "code", "literature_checks"))
from sim2e import poly_contact, random_rotations

def square(a=math.sqrt(math.pi)):
    h = a / 2
    return np.array([[-h, -h, 0], [h, -h, 0], [h, h, 0], [-h, h, 0]], float)

def triangle():
    R = math.sqrt(4 * math.pi / math.sqrt(3)) / math.sqrt(3)
    return np.array([[R * math.cos(t), R * math.sin(t), 0] for t in (math.pi / 2, math.pi / 2 + 2 * math.pi / 3, math.pi / 2 + 4 * math.pi / 3)])

def rz(t):
    c, s = math.cos(t), math.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])

def clip_area(P, Q):
    """area of the intersection of two convex counter-clockwise polygons in the plane (Sutherland-Hodgman)"""
    out = [tuple(p) for p in P]
    for k in range(len(Q)):
        a, b = Q[k], Q[(k + 1) % len(Q)]
        inp, out = out, []
        if not inp: break
        def side(p): return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
        for m in range(len(inp)):
            p, q = inp[m], inp[(m + 1) % len(inp)]
            sp, sq = side(p), side(q)
            if sp >= 0: out.append(p)
            if (sp > 0 and sq < 0) or (sp < 0 and sq > 0):
                f = sp / (sp - sq); out.append((p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1])))
    if len(out) < 3: return 0.0
    x = np.array(out); return 0.5 * abs(np.dot(x[:, 0], np.roll(x[:, 1], -1)) - np.dot(x[:, 1], np.roll(x[:, 0], -1)))

def to_xy(V, R, c):
    X = V @ R.T + c; return X[:, :2]

fails = 0
def check(name, got, want):
    global fails
    ok = (bool(got) == bool(want)); fails += (not ok)
    print(f"  {'OK  ' if ok else 'FAIL'} {name}: got {bool(got)} expect {bool(want)}")

I = np.eye(3)
print("coplanar plates")
for kind, V in (("square", square()), ("triangle", triangle())):
    rng = np.random.default_rng(3); bad = 0; both = 0
    for _ in range(4000):
        t = rng.uniform(0, 2 * math.pi); d = rng.uniform(-3, 3, 2)
        R = rz(t); dv = np.array([d[0], d[1], 0.0])
        got = poly_contact(I, dv, R, V)
        want = clip_area(to_xy(V, I, 0), to_xy(V, R, dv)) > 1e-9
        bad += (got != want); both += want
    print(f"  {kind}: 4000 random coplanar pairs, {both} overlapping, {bad} disagreements with the clipping oracle")
    fails += (bad > 0)
V = square()
check("coplanar overlap, same orientation", poly_contact(I, np.array([0.5, 0.3, 0]), I, V), True)
check("coplanar, separated", poly_contact(I, np.array([2.5, 0.0, 0]), I, V), False)
check("coplanar, edges touching only", poly_contact(I, np.array([math.sqrt(math.pi), 0.0, 0]), I, V), False)
check("parallel, plane gap 1e-3", poly_contact(I, np.array([0.1, 0.1, 1e-3]), I, V), False)
Rflip = np.diag([1.0, -1.0, -1.0])
check("antiparallel coplanar overlap", poly_contact(I, np.array([0.5, 0.3, 0]), Rflip, V), True)

print("vertex or edge in the other plate's plane: same answer from either side and with the roles swapped")
def tilted(axis, th):
    n = np.asarray(axis, float); n = n / np.linalg.norm(n)
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return I + math.sin(th) * K + (1 - math.cos(th)) * (K @ K)
for name, axis in (("edge", (1.0, 0.0, 0.0)), ("vertex", (1.0, 1.0, 0.0))):
    res = {}
    for flip in (1.0, -1.0):
        R = tilted(axis, flip * math.pi / 4)
        Xw = V @ R.T
        k = int(np.argmin(Xw[:, 2] * flip))          # the vertex that is to lie in the plane z = 0, at the origin
        dv = -Xw[k]                                  # all other vertices then lie on one side of the plane
        res[(flip, "ij")] = poly_contact(I, dv, R, V)
        res[(flip, "ji")] = poly_contact(R, -dv, I, V)
    vals = set(res.values())
    ok = len(vals) == 1
    fails += (not ok)
    print(f"  {'OK  ' if ok else 'FAIL'} {name} in plane, polygon on + side or - side, either role: {sorted(vals)}")
    if name == "edge":
        check("edge in plane lies inside the other plate: contact", vals.pop(), True) if ok else None

print("random non-parallel pairs against an independent chord test")
def indep(Ri, dv, Rj, V):
    ni, nj = Ri[:, 2], Rj[:, 2]; u = np.cross(ni, nj); u /= np.linalg.norm(u)
    def interval(R, c, pn, pc):
        X = V @ R.T + c; s = (X - pc) @ pn; ts = []
        for k in range(len(X)):
            k2 = (k + 1) % len(X)
            if (s[k] <= 0) != (s[k2] <= 0):
                f = s[k] / (s[k] - s[k2]); ts.append((X[k] + f * (X[k2] - X[k])) @ u)
        return (min(ts), max(ts)) if len(ts) >= 2 else None
    a = interval(Ri, np.zeros(3), nj, dv); b = interval(Rj, dv, ni, np.zeros(3))
    return a is not None and b is not None and max(a[0], b[0]) <= min(a[1], b[1])
rng = np.random.default_rng(5); bad = 0; hits = 0; n = 3000
for kind, V in (("square", square()), ("triangle", triangle())):
    Rs = random_rotations(2 * n, rng); bad = 0; hits = 0
    for k in range(n):
        dv = rng.uniform(-1.5, 1.5, 3)
        g = poly_contact(Rs[2 * k], dv, Rs[2 * k + 1], V); w = indep(Rs[2 * k], dv, Rs[2 * k + 1], V)
        bad += (g != w); hits += g
    print(f"  {kind}: {n} pairs, {hits} contacts, {bad} disagreements"); fails += (bad > 0)
print("FAILED" if fails else "ALL POLYGON CONTACT CHECKS PASSED")
sys.exit(1 if fails else 0)
