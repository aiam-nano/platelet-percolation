"""Self-tests required by INSTRUCTIONS.md steps 1 and 2. No reference thresholds."""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from numba import njit, prange

sys.path.insert(0, str(Path(__file__).resolve().parent))
import core

PASS = []
FAIL = []


def check(name, cond, detail=""):
    if cond:
        PASS.append(name)
        print(f"PASS  {name}" + (f"  {detail}" if detail else ""))
    else:
        FAIL.append(name)
        print(f"FAIL  {name}" + (f"  {detail}" if detail else ""))


def in_disc(x, p, n, R, tol=1e-7):
    d = x - p
    axial = abs(float(np.dot(d, n)))
    if axial > tol:
        return False
    return float(np.dot(d, d)) <= (R + tol) ** 2


def in_cyl(x, p, n, R, h, tol=1e-6):
    d = np.asarray(x, dtype=np.float64) - np.asarray(p, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    axial = float(np.dot(d, n))
    if abs(axial) > h + tol:
        return False
    rad2 = float(np.dot(d, d) - axial * axial)
    return rad2 <= (R + tol) ** 2


def ext_cyl(n, R, h, u):
    nu = abs(float(np.dot(n, u)))
    return h * nu + R * math.sqrt(max(0.0, 1.0 - nu * nu))


def in_ell(x, p, e1, e2, n, R, c, tol=1e-6):
    d = np.asarray(x) - np.asarray(p)
    a1 = float(np.dot(d, e1))
    a2 = float(np.dot(d, e2))
    an = float(np.dot(d, n))
    return (a1 / (R + tol)) ** 2 + (a2 / (R + tol)) ** 2 + (an / (c + tol)) ** 2 <= 1.0


def point_in_poly(pt, verts, nrm, tol=1e-6):
    if abs(float(np.dot(pt - verts[0], nrm))) > 1e-5:
        return False
    sign = 0.0
    nv = len(verts)
    for i in range(nv):
        j = (i + 1) % nv
        cr = np.cross(verts[j] - verts[i], pt - verts[i])
        s = float(np.dot(cr, nrm))
        if abs(s) <= tol:
            continue
        if sign == 0.0:
            sign = math.copysign(1.0, s)
        elif s * sign < 0.0:
            return False
    return True


def plate_verts(kind, p, e1, e2):
    if kind == 4:
        h = core.SQ_H
        loc = [(h, h), (-h, h), (-h, -h), (h, -h)]
    else:
        loc = [core.TRI_L0, core.TRI_L1, core.TRI_L2]
    return np.array([p + a * e1 + b * e2 for a, b in loc])


@njit(parallel=True, cache=True)
def _count_disc(pos, nrm):
    n = pos.shape[0]
    hits = np.zeros(n, np.int64)
    for i in prange(n):
        hits[i] = core.disc_hit(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
                                pos[i, 0], pos[i, 1], pos[i, 2],
                                nrm[i, 0], nrm[i, 1], nrm[i, 2], 1.0)
    return hits


@njit(parallel=True, cache=True)
def _count_cyl(pos, nrm, h):
    n = pos.shape[0]
    hits = np.zeros(n, np.int64)
    for i in prange(n):
        hits[i] = core.cylinder_hit(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, h,
                                    pos[i, 0], pos[i, 1], pos[i, 2],
                                    nrm[i, 0], nrm[i, 1], nrm[i, 2], 1.0, h)
    return hits


def random_orientations(rng, n):
    z = rng.uniform(-1.0, 1.0, size=n)
    phi = rng.uniform(0.0, 2.0 * math.pi, size=n)
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


def test_disc_witness_and_brute():
    rng = np.random.default_rng(1)
    n = 2000
    bad_witness = 0
    brute_disagree = 0
    sym_bad = 0
    for i in range(n):
        n1 = random_orientations(rng, 1)[0]
        n2 = random_orientations(rng, 1)[0]
        p2 = rng.normal(size=3)
        p2 *= rng.uniform(0.0, 2.2) / max(np.linalg.norm(p2), 1e-12)
        R1 = float(rng.uniform(0.5, 1.5))
        R2 = float(rng.uniform(0.5, 1.5))
        hit, wx, wy, wz = core.disc_analyze(0, 0, 0, n1[0], n1[1], n1[2], R1,
                                            p2[0], p2[1], p2[2], n2[0], n2[1], n2[2], R2)
        hit_b = core.disc_hit(p2[0], p2[1], p2[2], n2[0], n2[1], n2[2], R2,
                              0, 0, 0, n1[0], n1[1], n1[2], R1)
        if hit != hit_b:
            sym_bad += 1
        if hit:
            w = np.array([wx, wy, wz])
            if not in_disc(w, np.zeros(3), n1, R1) or not in_disc(w, p2, n2, R2):
                bad_witness += 1
        # Circumference samples: a hit there must agree.
        for host, centre, normal, radius, other_c, other_n, other_r in (
            (1, np.zeros(3), n1, R1, p2, n2, R2),
            (2, p2, n2, R2, np.zeros(3), n1, R1),
        ):
            e1, e2, normal = core.frame_from_normal(normal, 0.0)
            for k in range(64):
                ang = 2 * math.pi * k / 64
                x = centre + radius * (math.cos(ang) * e1 + math.sin(ang) * e2)
                if in_disc(x, other_c, other_n, other_r, tol=1e-8):
                    if not hit:
                        brute_disagree += 1
                    break
    check("disc witness in both discs", bad_witness == 0, f"bad={bad_witness}/{n}")
    check("disc symmetry", sym_bad == 0, f"bad={sym_bad}/{n}")
    check("disc circumference brute force", brute_disagree == 0, f"missed={brute_disagree}")


def test_disc_excluded():
    rng = np.random.default_rng(2)
    n = 400000
    pos = rng.uniform(-2.0, 2.0, size=(n, 3))
    nrm = random_orientations(rng, n)
    hits = int(_count_disc(pos, nrm).sum())
    vol = hits / n * 64.0
    # binomial SE
    p = hits / n
    se = 64.0 * math.sqrt(p * (1 - p) / n)
    target = math.pi ** 2
    check("disc excluded volume ~ pi^2", abs(vol - target) < 4 * se + 0.02,
          f"est={vol:.4f} se={se:.4f} target={target:.4f}")
    return vol, se


def test_cylinders():
    rng = np.random.default_rng(3)
    # Hand cases with a clear gap or overlap.
    cases = [
        # (p2, n2, h1, h2, expect) n1 = z, R = 1, p1 = 0
        ((0.5, 0, 0), (0, 0, 1), 0.2, 0.2, 1),
        ((2.5, 0, 0), (0, 0, 1), 0.2, 0.2, 0),
        ((0, 0, 0.3), (0, 0, 1), 0.2, 0.2, 1),
        ((0, 0, 0.6), (0, 0, 1), 0.2, 0.2, 0),
        ((0, 0, 0), (1, 0, 0), 0.5, 0.5, 1),
        ((0, 2.5, 0), (1, 0, 0), 0.5, 0.5, 0),
        ((0, 1.2, 0), (1, 0, 0), 0.5, 0.5, 1),
        ((3.0, 0, 0), (1, 0, 0), 0.5, 0.5, 0),
    ]
    hand_bad = 0
    for p2, n2, h1, h2, expect in cases:
        hit = core.cylinder_hit(0, 0, 0, 0, 0, 1, 1, h1, p2[0], p2[1], p2[2], n2[0], n2[1], n2[2], 1, h2)
        if hit != expect:
            hand_bad += 1
            print("  hand miss", p2, n2, h1, h2, "got", hit)
    check("cylinder hand cases", hand_bad == 0, f"bad={hand_bad}")

    n = 1500
    wit_bad = 0
    sep_bad = 0
    sym_bad = 0
    yes = 0
    for i in range(n):
        n1 = random_orientations(rng, 1)[0]
        n2 = random_orientations(rng, 1)[0]
        p2 = rng.normal(size=3)
        p2 *= rng.uniform(0.0, 2.4) / max(np.linalg.norm(p2), 1e-12)
        R1 = float(rng.uniform(0.4, 1.4))
        R2 = float(rng.uniform(0.4, 1.4))
        h1 = float(rng.uniform(0.02, 0.4))
        h2 = float(rng.uniform(0.02, 0.4))
        hit, wx, wy, wz, vx, vy, vz = core.cylinder_analyze(
            0, 0, 0, n1[0], n1[1], n1[2], R1, h1,
            p2[0], p2[1], p2[2], n2[0], n2[1], n2[2], R2, h2)
        hit2 = core.cylinder_hit(p2[0], p2[1], p2[2], n2[0], n2[1], n2[2], R2, h2,
                                 0, 0, 0, n1[0], n1[1], n1[2], R1, h1)
        if hit != hit2:
            sym_bad += 1
        if hit:
            yes += 1
            w = np.array([wx, wy, wz])
            if not in_cyl(w, np.zeros(3), n1, R1, h1) or not in_cyl(w, p2, n2, R2, h2):
                wit_bad += 1
        else:
            v = np.array([vx, vy, vz])
            dist = float(np.linalg.norm(v))
            if dist < 1e-12:
                sep_bad += 1
                continue
            u = v / dist
            gap = float(np.dot(-p2, u)) - ext_cyl(n1, R1, h1, u) - ext_cyl(n2, R2, h2, u)
            # (p1-p2)·u - ext1 - ext2 should equal |v|
            if abs(gap - dist) > 1e-5 or gap < -1e-6:
                sep_bad += 1
    check("cylinder witness", wit_bad == 0, f"bad={wit_bad} hits={yes}")
    check("cylinder separation identity", sep_bad == 0, f"bad={sep_bad}")
    check("cylinder symmetry", sym_bad == 0, f"bad={sym_bad}/{n}")


def test_cyl_disc_limit():
    rng = np.random.default_rng(4)
    n = 2000
    # Clearly overlapping or clearly separate discs, compared with a very thin cylinder.
    disagree = 0
    for i in range(n):
        n1 = random_orientations(rng, 1)[0]
        n2 = random_orientations(rng, 1)[0]
        p2 = rng.normal(size=3)
        p2 *= rng.uniform(0.0, 2.0) / max(np.linalg.norm(p2), 1e-12)
        hd, _, _, _ = core.disc_analyze(0, 0, 0, n1[0], n1[1], n1[2], 1,
                                        p2[0], p2[1], p2[2], n2[0], n2[1], n2[2], 1)
        hc = core.cylinder_hit(0, 0, 0, n1[0], n1[1], n1[2], 1, 0.0,
                               p2[0], p2[1], p2[2], n2[0], n2[1], n2[2], 1, 0.0)
        if hd != hc:
            # Ignore hairline cases: shift the centre by 1e-4 along p2 and see if the
            # disc answer is stable. Unstable answers sit on the boundary.
            if np.linalg.norm(p2) < 1e-8:
                continue
            direction = p2 / np.linalg.norm(p2)
            h_in, _, _, _ = core.disc_analyze(0, 0, 0, n1[0], n1[1], n1[2], 1,
                                             *(p2 - 1e-4 * direction), n2[0], n2[1], n2[2], 1)
            h_out, _, _, _ = core.disc_analyze(0, 0, 0, n1[0], n1[1], n1[2], 1,
                                              *(p2 + 1e-4 * direction), n2[0], n2[1], n2[2], 1)
            if h_in == h_out == hd:
                disagree += 1
    check("cylinder h=0 matches disc off the boundary", disagree == 0, f"disagree={disagree}/{n}")

    # Thickness tending to zero: fraction of disagreements against the disc test.
    fracs = []
    pos_n = 3000
    for h in (1e-2, 1e-3, 1e-4):
        bad = 0
        for i in range(pos_n):
            n1 = random_orientations(rng, 1)[0]
            n2 = random_orientations(rng, 1)[0]
            p2 = rng.uniform(-2, 2, size=3)
            hd = core.disc_hit(0, 0, 0, n1[0], n1[1], n1[2], 1, p2[0], p2[1], p2[2], n2[0], n2[1], n2[2], 1)
            hc = core.cylinder_hit(0, 0, 0, n1[0], n1[1], n1[2], 1, h, p2[0], p2[1], p2[2], n2[0], n2[1], n2[2], 1, h)
            if hd != hc:
                bad += 1
        fracs.append(bad / pos_n)
    decreasing = fracs[0] >= fracs[1] >= fracs[2]
    check("cylinder approaches disc as h -> 0", decreasing and fracs[-1] < 0.02,
          f"fracs={fracs}")
    return fracs


def test_cyl_excluded():
    # R=1, t=0.5, h=0.25. Analytic mean excluded volume from the specification.
    R = 1.0
    t = 0.5
    h = 0.25
    V = math.pi * R * R * t
    S = 2 * math.pi * R * t + 2 * math.pi * R * R
    M = math.pi * t + math.pi ** 2 * R
    target = 2 * V + S * M / (2 * math.pi)
    B = 2.0 * math.sqrt(R * R + h * h)
    rng = np.random.default_rng(5)
    n = 150000
    pos = rng.uniform(-B, B, size=(n, 3))
    nrm = random_orientations(rng, n)
    hits = int(_count_cyl(pos, nrm, h).sum())
    vol = hits / n * (2 * B) ** 3
    p = hits / n
    se = (2 * B) ** 3 * math.sqrt(max(p * (1 - p), 0.0) / n)
    check("cylinder excluded volume", abs(vol - target) < max(4 * se, 0.15),
          f"est={vol:.3f} se={se:.3f} target={target:.3f}")
    return vol, se, target


def test_plates():
    rng = np.random.default_rng(6)
    # Crossing squares through the origin.
    hit = core.plate_hit(4,
                         0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                         0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0)
    check("crossing squares intersect", hit == 1)
    # Same two squares translated well apart along x.
    far = core.SQ_SIDE + 0.2
    hit_far = core.plate_hit(4,
                             0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                             far, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0)
    check("separated squares miss", hit_far == 0)
    # Exactly coplanar overlapping squares are parallel: specified as no contact.
    cop = core.plate_hit(4,
                         0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                         0.2, 0.1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)
    check("coplanar squares reported as not intersecting", cop == 0)
    # Touching at a boundary point: centres separated by the side length along x,
    # both in xy, but that is coplanar. Use a perpendicular pair shifted to a vertex touch.
    # Square 1 in xy, square 2 in xz shifted by s along x: they share the point (s/2, 0, 0).
    s = core.SQ_SIDE
    touch = core.plate_hit(4,
                           0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                           s, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0)
    check("squares touching at a point", touch == 1, f"hit={touch}")

    wit_bad = 0
    sym_bad = 0
    n = 800
    for kind, nv in ((4, 4), (3, 3)):
        for i in range(n):
            n1 = random_orientations(rng, 1)[0]
            n2 = random_orientations(rng, 1)[0]
            e11, e21, n1 = core.frame_from_normal(n1, float(rng.uniform(0, 2 * math.pi)))
            e12, e22, n2 = core.frame_from_normal(n2, float(rng.uniform(0, 2 * math.pi)))
            p2 = rng.normal(size=3)
            p2 *= rng.uniform(0.0, 2.6) / max(np.linalg.norm(p2), 1e-12)
            hit, wx, wy, wz = core.plate_analyze(
                nv, 0, 0, 0, e11[0], e11[1], e11[2], e21[0], e21[1], e21[2], n1[0], n1[1], n1[2],
                p2[0], p2[1], p2[2], e12[0], e12[1], e12[2], e22[0], e22[1], e22[2], n2[0], n2[1], n2[2])
            hit2 = core.plate_hit(
                nv, p2[0], p2[1], p2[2], e12[0], e12[1], e12[2], e22[0], e22[1], e22[2], n2[0], n2[1], n2[2],
                0, 0, 0, e11[0], e11[1], e11[2], e21[0], e21[1], e21[2], n1[0], n1[1], n1[2])
            if hit != hit2:
                sym_bad += 1
            if hit:
                w = np.array([wx, wy, wz])
                v1 = plate_verts(nv, np.zeros(3), e11, e21)
                v2 = plate_verts(nv, p2, e12, e22)
                if not point_in_poly(w, v1, n1) or not point_in_poly(w, v2, n2):
                    wit_bad += 1
    check("plate witness", wit_bad == 0, f"bad={wit_bad}")
    check("plate symmetry", sym_bad == 0, f"bad={sym_bad}")


def _plate_exvol(kind_nv, B, n, seed):
    rng = np.random.default_rng(seed)
    pos = rng.uniform(-B, B, size=(n, 3))
    hits = 0
    # Fix plate 1 in the xy plane. Plate 2 gets a random frame.
    for i in range(n):
        n2 = random_orientations(rng, 1)[0]
        e1, e2, n2 = core.frame_from_normal(n2, float(rng.uniform(0, 2 * math.pi)))
        hits += core.plate_hit(kind_nv,
                               0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                               pos[i, 0], pos[i, 1], pos[i, 2],
                               e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], n2[0], n2[1], n2[2])
    vol = hits / n * (2 * B) ** 3
    p = hits / n
    se = (2 * B) ** 3 * math.sqrt(max(p * (1 - p), 0.0) / n)
    return vol, se


def test_plate_excluded():
    # A P / 2 for flat convex plates.
    A = math.pi
    P_sq = 4 * core.SQ_SIDE
    P_tri = 3 * core.TRI_SIDE
    t_sq = A * P_sq / 2
    t_tri = A * P_tri / 2
    # Bounding radii set the cube.
    v_sq, se_sq = _plate_exvol(4, core.RB_SQUARE * 2, 80000, 7)
    v_tri, se_tri = _plate_exvol(3, core.RB_TRIANGLE * 2, 80000, 8)
    check("square excluded volume", abs(v_sq - t_sq) < max(4 * se_sq, 0.2),
          f"est={v_sq:.3f} se={se_sq:.3f} target={t_sq:.3f}")
    check("triangle excluded volume", abs(v_tri - t_tri) < max(4 * se_tri, 0.25),
          f"est={v_tri:.3f} se={se_tri:.3f} target={t_tri:.3f}")
    return (v_sq, se_sq, t_sq), (v_tri, se_tri, t_tri)


def test_spheroids():
    rng = np.random.default_rng(9)
    n = 2000
    sphere_bad = 0
    near = 0
    wit_bad = 0
    sym_bad = 0
    for i in range(n):
        p2 = rng.normal(size=3)
        p2 *= rng.uniform(0.0, 3.2) / max(np.linalg.norm(p2), 1e-12)
        R1 = float(rng.uniform(0.5, 1.5))
        R2 = float(rng.uniform(0.5, 1.5))
        n1 = random_orientations(rng, 1)[0]
        n2 = random_orientations(rng, 1)[0]
        e11, e21, n1 = core.frame_from_normal(n1, float(rng.uniform(0, 6)))
        e12, e22, n2 = core.frame_from_normal(n2, float(rng.uniform(0, 6)))
        dist = float(np.linalg.norm(p2))
        expect = 1 if dist <= R1 + R2 else 0
        hit, wx, wy, wz, vx, vy, vz = core.spheroid_analyze(
            0, 0, 0, e11[0], e11[1], e11[2], e21[0], e21[1], e21[2], n1[0], n1[1], n1[2], R1, R1,
            p2[0], p2[1], p2[2], e12[0], e12[1], e12[2], e22[0], e22[1], e22[2], n2[0], n2[1], n2[2], R2, R2)
        if abs(dist - (R1 + R2)) <= 1e-6:
            near += 1
        elif hit != expect:
            sphere_bad += 1
        hit2 = core.spheroid_hit(
            p2[0], p2[1], p2[2], e12[0], e12[1], e12[2], e22[0], e22[1], e22[2], n2[0], n2[1], n2[2], R2, R2,
            0, 0, 0, e11[0], e11[1], e11[2], e21[0], e21[1], e21[2], n1[0], n1[1], n1[2], R1, R1)
        if hit != hit2:
            sym_bad += 1
        if hit:
            w = np.array([wx, wy, wz])
            if not in_ell(w, np.zeros(3), e11, e21, n1, R1, R1) or not in_ell(w, p2, e12, e22, n2, R2, R2):
                wit_bad += 1
    check("spheroid c=R matches spheres", sphere_bad == 0, f"bad={sphere_bad} near={near}")
    check("spheroid symmetry", sym_bad == 0, f"bad={sym_bad}")
    check("spheroid sphere witness", wit_bad == 0, f"bad={wit_bad}")

    # Flat limit. Count disagreements at several c.
    fracs = []
    m = 1500
    for c in (0.05, 0.01, 0.002):
        bad = 0
        for i in range(m):
            p2 = rng.uniform(-2, 2, size=3)
            n1 = random_orientations(rng, 1)[0]
            n2 = random_orientations(rng, 1)[0]
            e11, e21, n1 = core.frame_from_normal(n1, 0.3)
            e12, e22, n2 = core.frame_from_normal(n2, 1.1)
            hd = core.disc_hit(0, 0, 0, n1[0], n1[1], n1[2], 1.0,
                               p2[0], p2[1], p2[2], n2[0], n2[1], n2[2], 1.0)
            hs = core.spheroid_hit(
                0, 0, 0, e11[0], e11[1], e11[2], e21[0], e21[1], e21[2], n1[0], n1[1], n1[2], 1.0, c,
                p2[0], p2[1], p2[2], e12[0], e12[1], e12[2], e22[0], e22[1], e22[2], n2[0], n2[1], n2[2], 1.0, c)
            if hd != hs:
                bad += 1
        fracs.append(bad / m)
    check("spheroid approaches disc as c -> 0", fracs[0] >= fracs[1] >= fracs[2] and fracs[-1] < 0.05,
          f"fracs={fracs}")
    return fracs


def test_fisher():
    rng = np.random.default_rng(11)
    kappa = 6.312
    n = 200000
    acc = 0.0
    for i in range(n):
        c = core.fisher_cosine(rng.random(), kappa)
        acc += 0.5 * (3 * c * c - 1)
    mean = acc / n
    # Also the vector sampler about a fixed director, which is what the simulation uses.
    director = np.array([0.2, -0.4, 0.8])
    director /= np.linalg.norm(director)
    acc2 = 0.0
    m = 100000
    for i in range(m):
        v = core.fisher_normal(rng, director, kappa)
        c = float(np.dot(v, director))
        acc2 += 0.5 * (3 * c * c - 1)
    mean2 = acc2 / m
    check("Fisher <P2> about 0.60", abs(mean - 0.60) < 0.01 and abs(mean2 - 0.60) < 0.015,
          f"cosine={mean:.4f} vector={mean2:.4f}")
    return mean, mean2


def make_particles(rng, n, L, R=1.0):
    pos = rng.uniform(0.0, L, size=(n, 3))
    nrm = random_orientations(rng, n)
    e1 = np.zeros_like(pos)
    e2 = np.zeros_like(pos)
    rad = np.full(n, R)
    aux = np.zeros(n)
    rb = np.full(n, R)
    return pos, nrm, e1, e2, rad, aux, rb


def test_union_hand():
    L = 10.0
    xs = np.array([0.5, 2.0, 3.5, 5.0, 6.5, 8.0, 9.5])
    n = len(xs)
    pos = np.zeros((n, 3))
    pos[:, 0] = xs
    pos[:, 1] = 5.0
    pos[:, 2] = 5.0
    nrm = np.zeros((n, 3))
    nrm[:, 2] = 1.0
    e1 = np.zeros_like(pos)
    e2 = np.zeros_like(pos)
    e1[:, 0] = 1.0
    e2[:, 1] = 1.0
    rad = np.ones(n)
    aux = np.zeros(n)
    rb = np.ones(n)
    out = core.simulate(0, pos, nrm, e1, e2, rad, aux, rb, L, 4, 1)
    check("hand chain N_any and N_x", int(out[0]) == 7 and int(out[1]) == 7,
          f"out={out.tolist()}")
    check("hand chain no y or z wrap", int(out[2]) == -1 and int(out[3]) == -1, f"out={out.tolist()}")


def supercell_dirs(pos, nrm, rad, L, n_used):
    """Brute wrapping test: contact graph of the 3x3x3 supercell."""
    N = n_used
    parent = np.arange(N * 27, dtype=np.int64)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    def vid(i, ix, iy, iz):
        return i + N * ((ix + 1) + 3 * ((iy + 1) + 3 * (iz + 1)))

    offs = (-1, 0, 1)
    for i in range(N):
        for j in range(i + 1, N):
            for sx in offs:
                for sy in offs:
                    for sz in offs:
                        q = pos[j] + np.array([sx, sy, sz], dtype=np.float64) * L
                        lim = rad[i] + rad[j]
                        d = q - pos[i]
                        if float(np.dot(d, d)) > lim * lim:
                            continue
                        if not core.disc_hit(pos[i, 0], pos[i, 1], pos[i, 2],
                                             nrm[i, 0], nrm[i, 1], nrm[i, 2], rad[i],
                                             q[0], q[1], q[2],
                                             nrm[j, 0], nrm[j, 1], nrm[j, 2], rad[j]):
                            continue
                        for bx in offs:
                            for by in offs:
                                for bz in offs:
                                    cx, cy, cz = bx + sx, by + sy, bz + sz
                                    if -1 <= cx <= 1 and -1 <= cy <= 1 and -1 <= cz <= 1:
                                        union(vid(i, bx, by, bz), vid(j, cx, cy, cz))
    hx = hy = hz = False
    for i in range(N):
        r0 = find(vid(i, 0, 0, 0))
        for ax, ay, az in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
                           (1, 1, 0), (1, 0, 1), (0, 1, 1)):
            if find(vid(i, ax, ay, az)) == r0:
                if ax:
                    hx = True
                if ay:
                    hy = True
                if az:
                    hz = True
    return hx, hy, hz


def test_wrapping_and_cells():
    rng = np.random.default_rng(12)
    n_trials = 12
    bad = 0
    for t in range(n_trials):
        L = 12.0
        n = 36
        pos, nrm, e1, e2, rad, aux, rb = make_particles(rng, n, L, 1.0)
        # Compare grids. ncell=1 scans every particle, so it is the all-pairs search.
        a = core.simulate(0, pos, nrm, e1, e2, rad, aux, rb, L, 1, 1)
        b = core.simulate(0, pos, nrm, e1, e2, rad, aux, rb, L, 5, 1)
        c = core.simulate(0, pos, nrm, e1, e2, rad, aux, rb, L, 8, 1)
        if not (np.array_equal(a, b) and np.array_equal(a, c)):
            bad += 1
            print("  grid mismatch", a, b, c)
            continue
        # Brute supercell at N_any and just before, when a wrap was found.
        if a[0] > 0:
            n_at = int(a[0])
            hx, hy, hz = supercell_dirs(pos, nrm, rad, L, n_at)
            if (a[1] > 0) != hx or (a[2] > 0) != hy or (a[3] > 0) != hz:
                bad += 1
                print("  supercell mismatch at wrap", a, hx, hy, hz)
            if n_at > 1:
                hx0, hy0, hz0 = supercell_dirs(pos, nrm, rad, L, n_at - 1)
                if hx0 or hy0 or hz0:
                    bad += 1
                    print("  supercell wrapped early", a, hx0, hy0, hz0)
    check("cell list matches all-pairs and supercell", bad == 0, f"bad_trials={bad}/{n_trials}")


def main():
    print("triangle side", core.TRI_SIDE, "square side", core.SQ_SIDE)
    test_disc_witness_and_brute()
    vol = test_disc_excluded()
    test_cylinders()
    fr_c = test_cyl_disc_limit()
    cyl_ex = test_cyl_excluded()
    test_plates()
    plate_ex = test_plate_excluded()
    sph = test_spheroids()
    fisher = test_fisher()
    test_union_hand()
    test_wrapping_and_cells()
    print("---")
    print(f"{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        print("FAILED:", ", ".join(FAIL))
    results = {
        "disc_exvol": vol,
        "cyl_thin_fracs": fr_c,
        "cyl_exvol": cyl_ex,
        "plate_exvol": plate_ex,
        "sph_fracs": sph,
        "fisher": fisher,
        "fail": FAIL,
    }
    return results


if __name__ == "__main__":
    out = main()
    sys.exit(1 if FAIL else 0)
