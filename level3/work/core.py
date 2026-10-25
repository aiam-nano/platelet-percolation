"""Continuum percolation of penetrable platelets, cylinders and spheroids.

Geometry and the insertion loop. Choices that the specification leaves open
are recorded in NOTES.md.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

# Closed-contact tolerance in length units (radii are O(1)).
GJK_TOL = 1e-9
# Planes with sin^2(theta) below this are treated as parallel.
PARALLEL_SIN2 = 1e-24
# Vertex-on-plane tolerance for polygon clipping.
PLANE_EPS = 1e-12

SQ_SIDE = math.sqrt(math.pi)
SQ_H = 0.5 * SQ_SIDE
RB_SQUARE = SQ_SIDE / math.sqrt(2.0)
TRI_SIDE = math.sqrt(4.0 * math.pi / math.sqrt(3.0))
TRI_R = TRI_SIDE / math.sqrt(3.0)  # centroid to vertex
RB_TRIANGLE = TRI_R
_SQ3 = math.sqrt(3.0)
# Triangle vertices in (e1, e2), angles 90, 210, 330 degrees.
TRI_L0 = (TRI_R * 0.0, TRI_R * 1.0)
TRI_L1 = (TRI_R * (-0.5 * _SQ3), TRI_R * (-0.5))
TRI_L2 = (TRI_R * (0.5 * _SQ3), TRI_R * (-0.5))


@njit(cache=True, nogil=True)
def _cross(ax, ay, az, bx, by, bz):
    return ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx


@njit(cache=True, nogil=True)
def _dot(ax, ay, az, bx, by, bz):
    return ax * bx + ay * by + az * bz


@njit(cache=True, nogil=True)
def _triple(ax, ay, az, bx, by, bz, cx, cy, cz):
    return ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)


# ---------------------------------------------------------------------------
# Thin discs
# ---------------------------------------------------------------------------

@njit(cache=True, nogil=True)
def disc_analyze(p1x, p1y, p1z, n1x, n1y, n1z, R1,
                 p2x, p2y, p2z, n2x, n2y, n2z, R2):
    """Return (hit, wx, wy, wz). Witness is a point in both discs if hit."""
    cx, cy, cz = _cross(n1x, n1y, n1z, n2x, n2y, n2z)
    sin2 = cx * cx + cy * cy + cz * cz
    dx = p2x - p1x
    dy = p2y - p1y
    dz = p2z - p1z
    if sin2 < PARALLEL_SIN2:
        sep = abs(dx * n1x + dy * n1y + dz * n1z)
        # Nearly or exactly parallel. A positive plane separation means the
        # closed discs miss; a numerical zero is treated as coplanar.
        if sep > 1e-10 * (1.0 + R1 + R2):
            return 0, 0.0, 0.0, 0.0
        dist2 = dx * dx + dy * dy + dz * dz - sep * sep
        rad = R1 + R2
        if dist2 <= rad * rad:
            # Witness: a point on the segment joining the centres, clamped
            # into both discs. Centres lie in the common plane.
            if dist2 <= 1e-24:
                return 1, p1x, p1y, p1z
            dist = math.sqrt(dist2)
            # Move from p1 toward p2 by min(R1, dist) but stay within R2 of p2.
            # The overlap of two collinear intervals [0, dist] mapped to radii.
            # Point at distance t from p1 along the centre line, 0<=t<=dist,
            # needs t<=R1 and dist-t<=R2 i.e. t>=dist-R2.
            lo = 0.0
            if dist - R2 > lo:
                lo = dist - R2
            hi = dist
            if R1 < hi:
                hi = R1
            if lo <= hi:
                t = 0.5 * (lo + hi)
                inv = t / dist
                return 1, p1x + inv * dx, p1y + inv * dy, p1z + inv * dz
            return 1, p1x, p1y, p1z
        return 0, 0.0, 0.0, 0.0

    sinth = math.sqrt(sin2)
    inv_sin2 = 1.0 / sin2
    # Signed distance from p1 to plane 2, and from p2 to plane 1.
    dist1 = (p1x - p2x) * n2x + (p1y - p2y) * n2y + (p1z - p2z) * n2z
    dist2 = dx * n1x + dy * n1y + dz * n1z
    rho1 = abs(dist1) / sinth
    rho2 = abs(dist2) / sinth
    if rho1 > R1 or rho2 > R2:
        return 0, 0.0, 0.0, 0.0
    ux = cx / sinth
    uy = cy / sinth
    uz = cz / sinth
    n2dn1 = n2x * n1x + n2y * n1y + n2z * n1z
    # In-plane gradient of the other plane.
    v1x = n2x - n2dn1 * n1x
    v1y = n2y - n2dn1 * n1y
    v1z = n2z - n2dn1 * n1z
    f1x = p1x - dist1 * inv_sin2 * v1x
    f1y = p1y - dist1 * inv_sin2 * v1y
    f1z = p1z - dist1 * inv_sin2 * v1z
    v2x = n1x - n2dn1 * n2x
    v2y = n1y - n2dn1 * n2y
    v2z = n1z - n2dn1 * n2z
    f2x = p2x - dist2 * inv_sin2 * v2x
    f2y = p2y - dist2 * inv_sin2 * v2y
    f2z = p2z - dist2 * inv_sin2 * v2z
    s1 = math.sqrt(max(R1 * R1 - rho1 * rho1, 0.0))
    s2 = math.sqrt(max(R2 * R2 - rho2 * rho2, 0.0))
    t2 = (f2x - f1x) * ux + (f2y - f1y) * uy + (f2z - f1z) * uz
    a1 = -s1
    b1 = s1
    a2 = t2 - s2
    b2 = t2 + s2
    if a1 <= b2 and a2 <= b1:
        lo = a1 if a1 > a2 else a2
        hi = b1 if b1 < b2 else b2
        tm = 0.5 * (lo + hi)
        return 1, f1x + tm * ux, f1y + tm * uy, f1z + tm * uz
    return 0, 0.0, 0.0, 0.0


@njit(cache=True, nogil=True)
def disc_hit(p1x, p1y, p1z, n1x, n1y, n1z, R1,
             p2x, p2y, p2z, n2x, n2y, n2z, R2):
    hit, _, _, _ = disc_analyze(p1x, p1y, p1z, n1x, n1y, n1z, R1,
                                p2x, p2y, p2z, n2x, n2y, n2z, R2)
    return hit


# ---------------------------------------------------------------------------
# Convex polygons (square and triangular plates)
# ---------------------------------------------------------------------------

@njit(cache=True, nogil=True)
def _fill_square(px, py, pz, e1x, e1y, e1z, e2x, e2y, e2z, V):
    # Local corners (±h, ±h) in boundary order.
    h = SQ_H
    lx0 = h
    ly0 = h
    lx1 = -h
    ly1 = h
    lx2 = -h
    ly2 = -h
    lx3 = h
    ly3 = -h
    V[0, 0] = px + lx0 * e1x + ly0 * e2x
    V[0, 1] = py + lx0 * e1y + ly0 * e2y
    V[0, 2] = pz + lx0 * e1z + ly0 * e2z
    V[1, 0] = px + lx1 * e1x + ly1 * e2x
    V[1, 1] = py + lx1 * e1y + ly1 * e2y
    V[1, 2] = pz + lx1 * e1z + ly1 * e2z
    V[2, 0] = px + lx2 * e1x + ly2 * e2x
    V[2, 1] = py + lx2 * e1y + ly2 * e2y
    V[2, 2] = pz + lx2 * e1z + ly2 * e2z
    V[3, 0] = px + lx3 * e1x + ly3 * e2x
    V[3, 1] = py + lx3 * e1y + ly3 * e2y
    V[3, 2] = pz + lx3 * e1z + ly3 * e2z


@njit(cache=True, nogil=True)
def _fill_triangle(px, py, pz, e1x, e1y, e1z, e2x, e2y, e2z, V):
    lx0, ly0 = TRI_L0
    lx1, ly1 = TRI_L1
    lx2, ly2 = TRI_L2
    V[0, 0] = px + lx0 * e1x + ly0 * e2x
    V[0, 1] = py + lx0 * e1y + ly0 * e2y
    V[0, 2] = pz + lx0 * e1z + ly0 * e2z
    V[1, 0] = px + lx1 * e1x + ly1 * e2x
    V[1, 1] = py + lx1 * e1y + ly1 * e2y
    V[1, 2] = pz + lx1 * e1z + ly1 * e2z
    V[2, 0] = px + lx2 * e1x + ly2 * e2x
    V[2, 1] = py + lx2 * e1y + ly2 * e2y
    V[2, 2] = pz + lx2 * e1z + ly2 * e2z


@njit(cache=True, nogil=True)
def _chord(V, nv, nrx, nry, nrz, px, py, pz, ux, uy, uz):
    """Clip polygon against a plane. Return (nfound, tmin, tmax) along u."""
    d0 = (V[0, 0] - px) * nrx + (V[0, 1] - py) * nry + (V[0, 2] - pz) * nrz
    d1 = (V[1, 0] - px) * nrx + (V[1, 1] - py) * nry + (V[1, 2] - pz) * nrz
    d2 = (V[2, 0] - px) * nrx + (V[2, 1] - py) * nry + (V[2, 2] - pz) * nrz
    d3 = 0.0
    if nv == 4:
        d3 = (V[3, 0] - px) * nrx + (V[3, 1] - py) * nry + (V[3, 2] - pz) * nrz
    dists0 = d0
    dists1 = d1
    dists2 = d2
    dists3 = d3
    nfound = 0
    tmin = 0.0
    tmax = 0.0
    for i in range(nv):
        j = (i + 1) % nv
        if i == 0:
            di = dists0
        elif i == 1:
            di = dists1
        elif i == 2:
            di = dists2
        else:
            di = dists3
        if j == 0:
            dj = dists0
        elif j == 1:
            dj = dists1
        elif j == 2:
            dj = dists2
        else:
            dj = dists3
        if abs(di) <= PLANE_EPS:
            t = V[i, 0] * ux + V[i, 1] * uy + V[i, 2] * uz
            if nfound == 0:
                tmin = t
                tmax = t
            else:
                if t < tmin:
                    tmin = t
                if t > tmax:
                    tmax = t
            nfound += 1
        if (di > PLANE_EPS and dj < -PLANE_EPS) or (di < -PLANE_EPS and dj > PLANE_EPS):
            alpha = di / (di - dj)
            x = V[i, 0] + alpha * (V[j, 0] - V[i, 0])
            y = V[i, 1] + alpha * (V[j, 1] - V[i, 1])
            z = V[i, 2] + alpha * (V[j, 2] - V[i, 2])
            t = x * ux + y * uy + z * uz
            if nfound == 0:
                tmin = t
                tmax = t
            else:
                if t < tmin:
                    tmin = t
                if t > tmax:
                    tmax = t
            nfound += 1
    return nfound, tmin, tmax


@njit(cache=True, nogil=True)
def plate_analyze(nv, p1x, p1y, p1z, e11x, e11y, e11z, e21x, e21y, e21z, n1x, n1y, n1z,
                  p2x, p2y, p2z, e12x, e12y, e12z, e22x, e22y, e22z, n2x, n2y, n2z):
    """Zero-thickness convex plates. Parallel planes do not intersect."""
    cx, cy, cz = _cross(n1x, n1y, n1z, n2x, n2y, n2z)
    sin2 = cx * cx + cy * cy + cz * cz
    if sin2 < 1e-20:
        return 0, 0.0, 0.0, 0.0
    sinth = math.sqrt(sin2)
    ux = cx / sinth
    uy = cy / sinth
    uz = cz / sinth
    V1 = np.empty((4, 3))
    V2 = np.empty((4, 3))
    if nv == 4:
        _fill_square(p1x, p1y, p1z, e11x, e11y, e11z, e21x, e21y, e21z, V1)
        _fill_square(p2x, p2y, p2z, e12x, e12y, e12z, e22x, e22y, e22z, V2)
    else:
        _fill_triangle(p1x, p1y, p1z, e11x, e11y, e11z, e21x, e21y, e21z, V1)
        _fill_triangle(p2x, p2y, p2z, e12x, e12y, e12z, e22x, e22y, e22z, V2)
    nA, a1, b1 = _chord(V1, nv, n2x, n2y, n2z, p2x, p2y, p2z, ux, uy, uz)
    if nA == 0:
        return 0, 0.0, 0.0, 0.0
    nB, a2, b2 = _chord(V2, nv, n1x, n1y, n1z, p1x, p1y, p1z, ux, uy, uz)
    if nB == 0:
        return 0, 0.0, 0.0, 0.0
    if a1 <= b2 and a2 <= b1:
        lo = a1 if a1 > a2 else a2
        hi = b1 if b1 < b2 else b2
        tm = 0.5 * (lo + hi)
        # A point on the line: reconstruct from plane intersection using p1's foot.
        # Use the projection origin at 0, so the point is tm * u plus the component
        # of any chord point orthogonal to u. Recover from polygon 1's chord midpoint
        # by solving: point = q + tm_rel * u, but tm is absolute dot(point, u).
        # point = (component of a known point on the line). Take a point on chord A
        # with parameter a1, then step (tm - a1) along u.
        # Rebuild one endpoint at parameter a1 by clipping again is wasteful.
        # The line is not through the origin, so point = P + ((tm - dot(P,u)) * u)
        # where P is any point on the line. Use the foot of p1 if it lies on the
        # line, which it does only if p1 is on plane 2. Instead compute P as the
        # intersection point closest to p1, same formula as the disc foot.
        dist1 = (p1x - p2x) * n2x + (p1y - p2y) * n2y + (p1z - p2z) * n2z
        n2dn1 = n2x * n1x + n2y * n1y + n2z * n1z
        v1x = n2x - n2dn1 * n1x
        v1y = n2y - n2dn1 * n1y
        v1z = n2z - n2dn1 * n1z
        inv_sin2 = 1.0 / sin2
        f1x = p1x - dist1 * inv_sin2 * v1x
        f1y = p1y - dist1 * inv_sin2 * v1y
        f1z = p1z - dist1 * inv_sin2 * v1z
        tfoot = f1x * ux + f1y * uy + f1z * uz
        step = tm - tfoot
        return 1, f1x + step * ux, f1y + step * uy, f1z + step * uz
    return 0, 0.0, 0.0, 0.0


@njit(cache=True, nogil=True)
def plate_hit(nv, p1x, p1y, p1z, e11x, e11y, e11z, e21x, e21y, e21z, n1x, n1y, n1z,
              p2x, p2y, p2z, e12x, e12y, e12z, e22x, e22y, e22z, n2x, n2y, n2z):
    hit, _, _, _ = plate_analyze(
        nv, p1x, p1y, p1z, e11x, e11y, e11z, e21x, e21y, e21z, n1x, n1y, n1z,
        p2x, p2y, p2z, e12x, e12y, e12z, e22x, e22y, e22z, n2x, n2y, n2z)
    return hit


# ---------------------------------------------------------------------------
# Support functions and GJK
# ---------------------------------------------------------------------------

@njit(cache=True, nogil=True)
def _support_cyl(px, py, pz, nx, ny, nz, R, h, dx, dy, dz):
    vn = dx * nx + dy * ny + dz * nz
    if vn >= 0.0:
        ax = h * nx
        ay = h * ny
        az = h * nz
    else:
        ax = -h * nx
        ay = -h * ny
        az = -h * nz
    vx = dx - vn * nx
    vy = dy - vn * ny
    vz = dz - vn * nz
    vlen = math.sqrt(vx * vx + vy * vy + vz * vz)
    if vlen > 1e-16:
        s = R / vlen
        rx = vx * s
        ry = vy * s
        rz = vz * s
    else:
        rx = 0.0
        ry = 0.0
        rz = 0.0
    return px + ax + rx, py + ay + ry, pz + az + rz


@njit(cache=True, nogil=True)
def _support_ell(px, py, pz, e1x, e1y, e1z, e2x, e2y, e2z, nx, ny, nz, R, c, dx, dy, dz):
    a1 = dx * e1x + dy * e1y + dz * e1z
    a2 = dx * e2x + dy * e2y + dz * e2z
    an = dx * nx + dy * ny + dz * nz
    mnorm = math.sqrt((R * a1) * (R * a1) + (R * a2) * (R * a2) + (c * an) * (c * an))
    if mnorm < 1e-30:
        return px, py, pz
    inv = 1.0 / mnorm
    return (px + (R * R * a1 * inv) * e1x + (R * R * a2 * inv) * e2x + (c * c * an * inv) * nx,
            py + (R * R * a1 * inv) * e1y + (R * R * a2 * inv) * e2y + (c * c * an * inv) * ny,
            pz + (R * R * a1 * inv) * e1z + (R * R * a2 * inv) * e2z + (c * c * an * inv) * nz)


@njit(cache=True, nogil=True)
def _closest_segment(ax, ay, az, bx, by, bz):
    abx = bx - ax
    aby = by - ay
    abz = bz - az
    den = abx * abx + aby * aby + abz * abz
    if den < 1e-30:
        return ax, ay, az, 1.0, 0.0
    t = -(ax * abx + ay * aby + az * abz) / den
    if t <= 0.0:
        return ax, ay, az, 1.0, 0.0
    if t >= 1.0:
        return bx, by, bz, 0.0, 1.0
    return ax + t * abx, ay + t * aby, az + t * abz, 1.0 - t, t


@njit(cache=True, nogil=True)
def _closest_triangle(ax, ay, az, bx, by, bz, cx, cy, cz):
    """Closest point of triangle ABC to the origin, plus barycentric weights of A,B,C."""
    abx = bx - ax
    aby = by - ay
    abz = bz - az
    acx = cx - ax
    acy = cy - ay
    acz = cz - az
    # d1 = AB·(O-A) = -AB·A
    d1 = -(abx * ax + aby * ay + abz * az)
    d2 = -(acx * ax + acy * ay + acz * az)
    if d1 <= 0.0 and d2 <= 0.0:
        return ax, ay, az, 1.0, 0.0, 0.0
    d3 = -(abx * bx + aby * by + abz * bz)
    d4 = -(acx * bx + acy * by + acz * bz)
    if d3 >= 0.0 and d4 <= d3:
        return bx, by, bz, 0.0, 1.0, 0.0
    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        den = d1 - d3
        v = d1 / den if den != 0.0 else 0.0
        return ax + v * abx, ay + v * aby, az + v * abz, 1.0 - v, v, 0.0
    d5 = -(abx * cx + aby * cy + abz * cz)
    d6 = -(acx * cx + acy * cy + acz * cz)
    if d6 >= 0.0 and d5 <= d6:
        return cx, cy, cz, 0.0, 0.0, 1.0
    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        den = d2 - d6
        w = d2 / den if den != 0.0 else 0.0
        return ax + w * acx, ay + w * acy, az + w * acz, 1.0 - w, 0.0, w
    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        den = (d4 - d3) + (d5 - d6)
        w = (d4 - d3) / den if den != 0.0 else 0.0
        return (bx + w * (cx - bx), by + w * (cy - by), bz + w * (cz - bz),
                0.0, 1.0 - w, w)
    denom = va + vb + vc
    if denom <= 1e-30:
        # Degenerate face: pick the closest vertex.
        da = ax * ax + ay * ay + az * az
        db = bx * bx + by * by + bz * bz
        dc = cx * cx + cy * cy + cz * cz
        if da <= db and da <= dc:
            return ax, ay, az, 1.0, 0.0, 0.0
        if db <= dc:
            return bx, by, bz, 0.0, 1.0, 0.0
        return cx, cy, cz, 0.0, 0.0, 1.0
    v = vb / denom
    w = vc / denom
    u = 1.0 - v - w
    return ax + v * abx + w * acx, ay + v * aby + w * acy, az + v * abz + w * acz, u, v, w


@njit(cache=True, nogil=True)
def _tet_bary(ax, ay, az, bx, by, bz, cx, cy, cz, ddx, ddy, ddz, vx, vy, vz):
    """Barycentric weights of point v in tet ABCD. Returns (ok, wa, wb, wc, wd)."""
    c0x = bx - ax
    c0y = by - ay
    c0z = bz - az
    c1x = cx - ax
    c1y = cy - ay
    c1z = cz - az
    c2x = ddx - ax
    c2y = ddy - ay
    c2z = ddz - az
    rhsx = vx - ax
    rhsy = vy - ay
    rhsz = vz - az
    det = _triple(c0x, c0y, c0z, c1x, c1y, c1z, c2x, c2y, c2z)
    if abs(det) < 1e-18:
        return 0, 0.0, 0.0, 0.0, 0.0
    wb = _triple(rhsx, rhsy, rhsz, c1x, c1y, c1z, c2x, c2y, c2z) / det
    wc = _triple(c0x, c0y, c0z, rhsx, rhsy, rhsz, c2x, c2y, c2z) / det
    wd = _triple(c0x, c0y, c0z, c1x, c1y, c1z, rhsx, rhsy, rhsz) / det
    wa = 1.0 - wb - wc - wd
    return 1, wa, wb, wc, wd


@njit(cache=True, nogil=True)
def _gjk(kind,
         p1x, p1y, p1z, n1x, n1y, n1z, e11x, e11y, e11z, e21x, e21y, e21z, R1, a1,
         p2x, p2y, p2z, n2x, n2y, n2z, e12x, e12y, e12z, e22x, e22y, e22z, R2, a2):
    """GJK on cylinder (kind=0, a=half-height) or ellipsoid (kind=1, a=polar semi-axis).

    Returns hit, witness xyz, closest-point-of-Minkowski xyz.
    """
    # Zeros, not empty: a zero barycentric weight must not touch uninitialised memory.
    pts = np.zeros((4, 3))
    pa = np.zeros((4, 3))
    pb = np.zeros((4, 3))
    dx = p1x - p2x
    dy = p1y - p2y
    dz = p1z - p2z
    if dx * dx + dy * dy + dz * dz < 1e-24:
        dx, dy, dz = 1.0, 0.0, 0.0

    def support(sx, sy, sz):
        if kind == 0:
            ax, ay, az = _support_cyl(p1x, p1y, p1z, n1x, n1y, n1z, R1, a1, sx, sy, sz)
            bx, by, bz = _support_cyl(p2x, p2y, p2z, n2x, n2y, n2z, R2, a2, -sx, -sy, -sz)
        else:
            ax, ay, az = _support_ell(p1x, p1y, p1z, e11x, e11y, e11z, e21x, e21y, e21z,
                                      n1x, n1y, n1z, R1, a1, sx, sy, sz)
            bx, by, bz = _support_ell(p2x, p2y, p2z, e12x, e12y, e12z, e22x, e22y, e22z,
                                      n2x, n2y, n2z, R2, a2, -sx, -sy, -sz)
        return ax - bx, ay - by, az - bz, ax, ay, az, bx, by, bz

    wx, wy, wz, pax, pay, paz, pbx, pby, pbz = support(dx, dy, dz)
    pts[0, 0] = wx
    pts[0, 1] = wy
    pts[0, 2] = wz
    pa[0, 0] = pax
    pa[0, 1] = pay
    pa[0, 2] = paz
    pb[0, 0] = pbx
    pb[0, 1] = pby
    pb[0, 2] = pbz
    n = 1
    tol = GJK_TOL
    tol2 = tol * tol

    for _it in range(64):
        # Reduce to the feature closest to the origin and get barycentric weights
        # on the UNREDUCED indexing, then compact.
        if n == 1:
            vx = pts[0, 0]
            vy = pts[0, 1]
            vz = pts[0, 2]
            w0, w1, w2, w3 = 1.0, 0.0, 0.0, 0.0
            nred = 1
        elif n == 2:
            vx, vy, vz, wa, wb = _closest_segment(
                pts[0, 0], pts[0, 1], pts[0, 2], pts[1, 0], pts[1, 1], pts[1, 2])
            w0, w1, w2, w3 = wa, wb, 0.0, 0.0
            # Compact to the vertices with nonzero weight.
            if wa > 1e-14 and wb > 1e-14:
                nred = 2
            elif wb > 1e-14:
                pts[0, 0] = pts[1, 0]
                pts[0, 1] = pts[1, 1]
                pts[0, 2] = pts[1, 2]
                pa[0, 0] = pa[1, 0]
                pa[0, 1] = pa[1, 1]
                pa[0, 2] = pa[1, 2]
                pb[0, 0] = pb[1, 0]
                pb[0, 1] = pb[1, 1]
                pb[0, 2] = pb[1, 2]
                w0, w1 = 1.0, 0.0
                nred = 1
            else:
                w0, w1 = 1.0, 0.0
                nred = 1
        elif n == 3:
            vx, vy, vz, wa, wb, wc = _closest_triangle(
                pts[0, 0], pts[0, 1], pts[0, 2],
                pts[1, 0], pts[1, 1], pts[1, 2],
                pts[2, 0], pts[2, 1], pts[2, 2])
            w0, w1, w2, w3 = wa, wb, wc, 0.0
            mask0 = 1 if wa > 1e-14 else 0
            mask1 = 1 if wb > 1e-14 else 0
            mask2 = 1 if wc > 1e-14 else 0
            nred = mask0 + mask1 + mask2
            if nred == 0:
                nred = 1
                mask0 = 1
                w0 = 1.0
            if nred < 3:
                tmpx0, tmpy0, tmpz0 = pts[0, 0], pts[0, 1], pts[0, 2]
                tmpx1, tmpy1, tmpz1 = pts[1, 0], pts[1, 1], pts[1, 2]
                tmpx2, tmpy2, tmpz2 = pts[2, 0], pts[2, 1], pts[2, 2]
                pax0, pay0, paz0 = pa[0, 0], pa[0, 1], pa[0, 2]
                pax1, pay1, paz1 = pa[1, 0], pa[1, 1], pa[1, 2]
                pax2, pay2, paz2 = pa[2, 0], pa[2, 1], pa[2, 2]
                pbx0, pby0, pbz0 = pb[0, 0], pb[0, 1], pb[0, 2]
                pbx1, pby1, pbz1 = pb[1, 0], pb[1, 1], pb[1, 2]
                pbx2, pby2, pbz2 = pb[2, 0], pb[2, 1], pb[2, 2]
                ww = (w0, w1, w2)
                xs = (tmpx0, tmpx1, tmpx2)
                ys = (tmpy0, tmpy1, tmpy2)
                zs = (tmpz0, tmpz1, tmpz2)
                pas = ((pax0, pay0, paz0), (pax1, pay1, paz1), (pax2, pay2, paz2))
                pbs = ((pbx0, pby0, pbz0), (pbx1, pby1, pbz1), (pbx2, pby2, pbz2))
                masks = (mask0, mask1, mask2)
                k = 0
                nw0 = 0.0
                nw1 = 0.0
                nw2 = 0.0
                for t in range(3):
                    if masks[t]:
                        pts[k, 0] = xs[t]
                        pts[k, 1] = ys[t]
                        pts[k, 2] = zs[t]
                        pa[k, 0] = pas[t][0]
                        pa[k, 1] = pas[t][1]
                        pa[k, 2] = pas[t][2]
                        pb[k, 0] = pbs[t][0]
                        pb[k, 1] = pbs[t][1]
                        pb[k, 2] = pbs[t][2]
                        if k == 0:
                            nw0 = ww[t]
                        elif k == 1:
                            nw1 = ww[t]
                        else:
                            nw2 = ww[t]
                        k += 1
                w0, w1, w2 = nw0, nw1, nw2
                nred = k
        else:
            ok, wa, wb, wc, wd = _tet_bary(
                pts[0, 0], pts[0, 1], pts[0, 2],
                pts[1, 0], pts[1, 1], pts[1, 2],
                pts[2, 0], pts[2, 1], pts[2, 2],
                pts[3, 0], pts[3, 1], pts[3, 2],
                0.0, 0.0, 0.0)
            inside = (ok == 1 and wa >= -1e-8 and wb >= -1e-8 and wc >= -1e-8 and wd >= -1e-8)
            if inside:
                vx, vy, vz = 0.0, 0.0, 0.0
                w0, w1, w2, w3 = wa, wb, wc, wd
                nred = 4
            else:
                best = 1e300
                vx = pts[0, 0]
                vy = pts[0, 1]
                vz = pts[0, 2]
                best_face = 0
                bw0 = 1.0
                bw1 = 0.0
                bw2 = 0.0
                faces0 = (0, 0, 0, 1)
                faces1 = (1, 1, 2, 2)
                faces2 = (2, 3, 3, 3)
                for f in range(4):
                    ia = faces0[f]
                    ib = faces1[f]
                    ic = faces2[f]
                    cx1, cy1, cz1, ua, ub, uc = _closest_triangle(
                        pts[ia, 0], pts[ia, 1], pts[ia, 2],
                        pts[ib, 0], pts[ib, 1], pts[ib, 2],
                        pts[ic, 0], pts[ic, 1], pts[ic, 2])
                    d2 = cx1 * cx1 + cy1 * cy1 + cz1 * cz1
                    if d2 < best:
                        best = d2
                        vx, vy, vz = cx1, cy1, cz1
                        best_face = f
                        bw0, bw1, bw2 = ua, ub, uc
                ia = faces0[best_face]
                ib = faces1[best_face]
                ic = faces2[best_face]
                # Compact the face, dropping near-zero weights.
                ids = (ia, ib, ic)
                bws = (bw0, bw1, bw2)
                k = 0
                nw = (0.0, 0.0, 0.0)
                # manual copy
                sel = np.empty(3, np.int64)
                selw = np.empty(3)
                for t in range(3):
                    if bws[t] > 1e-14:
                        sel[k] = ids[t]
                        selw[k] = bws[t]
                        k += 1
                if k == 0:
                    sel[0] = ia
                    selw[0] = 1.0
                    k = 1
                # stash
                tx = np.empty(3)
                ty = np.empty(3)
                tz = np.empty(3)
                tax = np.empty(3)
                tay = np.empty(3)
                taz = np.empty(3)
                tbx = np.empty(3)
                tby = np.empty(3)
                tbz = np.empty(3)
                for t in range(k):
                    i = sel[t]
                    tx[t] = pts[i, 0]
                    ty[t] = pts[i, 1]
                    tz[t] = pts[i, 2]
                    tax[t] = pa[i, 0]
                    tay[t] = pa[i, 1]
                    taz[t] = pa[i, 2]
                    tbx[t] = pb[i, 0]
                    tby[t] = pb[i, 1]
                    tbz[t] = pb[i, 2]
                for t in range(k):
                    pts[t, 0] = tx[t]
                    pts[t, 1] = ty[t]
                    pts[t, 2] = tz[t]
                    pa[t, 0] = tax[t]
                    pa[t, 1] = tay[t]
                    pa[t, 2] = taz[t]
                    pb[t, 0] = tbx[t]
                    pb[t, 1] = tby[t]
                    pb[t, 2] = tbz[t]
                w0 = selw[0]
                w1 = selw[1] if k > 1 else 0.0
                w2 = selw[2] if k > 2 else 0.0
                w3 = 0.0
                nred = k
        n = nred
        dist2 = vx * vx + vy * vy + vz * vz
        if dist2 <= tol2:
            wax = w0 * pa[0, 0] + w1 * pa[1, 0] + w2 * pa[2, 0] + w3 * pa[3, 0]
            way = w0 * pa[0, 1] + w1 * pa[1, 1] + w2 * pa[2, 1] + w3 * pa[3, 1]
            waz = w0 * pa[0, 2] + w1 * pa[1, 2] + w2 * pa[2, 2] + w3 * pa[3, 2]
            return 1, wax, way, waz, vx, vy, vz
        # Search toward the origin.
        sx = -vx
        sy = -vy
        sz = -vz
        mx, my, mz, pax, pay, paz, pbx, pby, pbz = support(sx, sy, sz)
        # Improvement of the lower bound on distance.
        improve = dist2 - (mx * vx + my * vy + mz * vz)
        if improve <= tol * math.sqrt(dist2):
            return 0, 0.0, 0.0, 0.0, vx, vy, vz
        # Reject a duplicate support point.
        dup = False
        for i in range(n):
            ddx = mx - pts[i, 0]
            ddy = my - pts[i, 1]
            ddz = mz - pts[i, 2]
            if ddx * ddx + ddy * ddy + ddz * ddz < 1e-24:
                dup = True
                break
        if dup:
            return 0, 0.0, 0.0, 0.0, vx, vy, vz
        pts[n, 0] = mx
        pts[n, 1] = my
        pts[n, 2] = mz
        pa[n, 0] = pax
        pa[n, 1] = pay
        pa[n, 2] = paz
        pb[n, 0] = pbx
        pb[n, 1] = pby
        pb[n, 2] = pbz
        n += 1
    # Iteration cap: treat a tiny residual as contact.
    if dist2 <= tol2:
        wax = w0 * pa[0, 0] + w1 * pa[1, 0] + w2 * pa[2, 0] + w3 * pa[3, 0]
        way = w0 * pa[0, 1] + w1 * pa[1, 1] + w2 * pa[2, 1] + w3 * pa[3, 1]
        waz = w0 * pa[0, 2] + w1 * pa[1, 2] + w2 * pa[2, 2] + w3 * pa[3, 2]
        return 1, wax, way, waz, vx, vy, vz
    return 0, 0.0, 0.0, 0.0, vx, vy, vz


@njit(cache=True, nogil=True)
def cylinder_analyze(p1x, p1y, p1z, n1x, n1y, n1z, R1, h1,
                     p2x, p2y, p2z, n2x, n2y, n2z, R2, h2):
    # Bounding-sphere reject. The cylinder fits in a ball of radius sqrt(R^2+h^2).
    dx = p2x - p1x
    dy = p2y - p1y
    dz = p2z - p1z
    rb1 = math.sqrt(R1 * R1 + h1 * h1)
    rb2 = math.sqrt(R2 * R2 + h2 * h2)
    if dx * dx + dy * dy + dz * dz > (rb1 + rb2) * (rb1 + rb2):
        return 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    return _gjk(0,
                p1x, p1y, p1z, n1x, n1y, n1z, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, R1, h1,
                p2x, p2y, p2z, n2x, n2y, n2z, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, R2, h2)


@njit(cache=True, nogil=True)
def cylinder_hit(p1x, p1y, p1z, n1x, n1y, n1z, R1, h1,
                 p2x, p2y, p2z, n2x, n2y, n2z, R2, h2):
    hit, _, _, _, _, _, _ = cylinder_analyze(
        p1x, p1y, p1z, n1x, n1y, n1z, R1, h1,
        p2x, p2y, p2z, n2x, n2y, n2z, R2, h2)
    return hit


@njit(cache=True, nogil=True)
def spheroid_analyze(p1x, p1y, p1z, e11x, e11y, e11z, e21x, e21y, e21z, n1x, n1y, n1z, R1, c1,
                     p2x, p2y, p2z, e12x, e12y, e12z, e22x, e22y, e22z, n2x, n2y, n2z, R2, c2):
    dx = p2x - p1x
    dy = p2y - p1y
    dz = p2z - p1z
    rb1 = R1 if R1 > c1 else c1
    rb2 = R2 if R2 > c2 else c2
    if dx * dx + dy * dy + dz * dz > (rb1 + rb2) * (rb1 + rb2):
        return 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    return _gjk(1,
                p1x, p1y, p1z, n1x, n1y, n1z, e11x, e11y, e11z, e21x, e21y, e21z, R1, c1,
                p2x, p2y, p2z, n2x, n2y, n2z, e12x, e12y, e12z, e22x, e22y, e22z, R2, c2)


@njit(cache=True, nogil=True)
def spheroid_hit(p1x, p1y, p1z, e11x, e11y, e11z, e21x, e21y, e21z, n1x, n1y, n1z, R1, c1,
                 p2x, p2y, p2z, e12x, e12y, e12z, e22x, e22y, e22z, n2x, n2y, n2z, R2, c2):
    hit, _, _, _, _, _, _ = spheroid_analyze(
        p1x, p1y, p1z, e11x, e11y, e11z, e21x, e21y, e21z, n1x, n1y, n1z, R1, c1,
        p2x, p2y, p2z, e12x, e12y, e12z, e22x, e22y, e22z, n2x, n2y, n2z, R2, c2)
    return hit


# ---------------------------------------------------------------------------
# Insertion, cell list, periodic union-find
# ---------------------------------------------------------------------------

@njit(cache=True, nogil=True)
def _find(x, parent, shift):
    """Iterative find with offset compression.

    shift[i] is the lattice image of i that sits in the root's cell 0.
    """
    stack = np.empty(256, np.int64)
    sp = 0
    i = x
    while parent[i] != i:
        stack[sp] = i
        sp += 1
        i = parent[i]
        if sp >= 256:
            break
    root = i
    while sp > 0:
        sp -= 1
        node = stack[sp]
        p = parent[node]
        shift[node, 0] += shift[p, 0]
        shift[node, 1] += shift[p, 1]
        shift[node, 2] += shift[p, 2]
        parent[node] = root
    return root


@njit(cache=True, nogil=True)
def _union(x, y, s0, s1, s2, parent, shift, rank):
    """Contact between x at image 0 and y at image s. Return winding xyz.

    A nonzero winding means the contact closes a wrapping cycle. That edge is
    not unioned, so the offset labelling stays single-valued.
    """
    rx = _find(x, parent, shift)
    sx0 = shift[x, 0]
    sx1 = shift[x, 1]
    sx2 = shift[x, 2]
    ry = _find(y, parent, shift)
    sy0 = shift[y, 0]
    sy1 = shift[y, 1]
    sy2 = shift[y, 2]
    if rx == ry:
        return s0 + sx0 - sy0, s1 + sx1 - sy1, s2 + sx2 - sy2
    # Attach the lower rank root under the higher one.
    if rank[rx] < rank[ry]:
        parent[rx] = ry
        shift[rx, 0] = -s0 - sx0 + sy0
        shift[rx, 1] = -s1 - sx1 + sy1
        shift[rx, 2] = -s2 - sx2 + sy2
    else:
        parent[ry] = rx
        shift[ry, 0] = s0 + sx0 - sy0
        shift[ry, 1] = s1 + sx1 - sy1
        shift[ry, 2] = s2 + sx2 - sy2
        if rank[rx] == rank[ry]:
            rank[rx] += 1
    return 0, 0, 0


@njit(cache=True, nogil=True)
def _min_image_shift(pix, piy, piz, pjx, pjy, pjz, L):
    s0 = int(-math.floor((pjx - pix) / L + 0.5))
    s1 = int(-math.floor((pjy - piy) / L + 0.5))
    s2 = int(-math.floor((pjz - piz) / L + 0.5))
    return s0, s1, s2


@njit(cache=True, nogil=True)
def _cell_of(p, L, ncell):
    u = p / L
    if u < 0.0:
        u = 0.0
    if u >= 1.0:
        u = 0.0
    c = int(u * ncell)
    if c >= ncell:
        c = ncell - 1
    if c < 0:
        c = 0
    return c


@njit(cache=True, nogil=True)
def _contact(kind, i, j, qjx, qjy, qjz,
             pos, nrm, e1, e2, rad, aux):
    """Contact of particle i at its stored centre with particle j translated to qj."""
    if kind == 0:
        return disc_hit(pos[i, 0], pos[i, 1], pos[i, 2], nrm[i, 0], nrm[i, 1], nrm[i, 2], rad[i],
                        qjx, qjy, qjz, nrm[j, 0], nrm[j, 1], nrm[j, 2], rad[j])
    if kind == 1:
        return cylinder_hit(pos[i, 0], pos[i, 1], pos[i, 2], nrm[i, 0], nrm[i, 1], nrm[i, 2], rad[i], aux[i],
                            qjx, qjy, qjz, nrm[j, 0], nrm[j, 1], nrm[j, 2], rad[j], aux[j])
    if kind == 2 or kind == 3:
        nv = 4 if kind == 2 else 3
        return plate_hit(nv,
                         pos[i, 0], pos[i, 1], pos[i, 2],
                         e1[i, 0], e1[i, 1], e1[i, 2], e2[i, 0], e2[i, 1], e2[i, 2],
                         nrm[i, 0], nrm[i, 1], nrm[i, 2],
                         qjx, qjy, qjz,
                         e1[j, 0], e1[j, 1], e1[j, 2], e2[j, 0], e2[j, 1], e2[j, 2],
                         nrm[j, 0], nrm[j, 1], nrm[j, 2])
    # kind == 4 spheroid; aux is the polar semi-axis c
    return spheroid_hit(pos[i, 0], pos[i, 1], pos[i, 2],
                        e1[i, 0], e1[i, 1], e1[i, 2], e2[i, 0], e2[i, 1], e2[i, 2],
                        nrm[i, 0], nrm[i, 1], nrm[i, 2], rad[i], aux[i],
                        qjx, qjy, qjz,
                        e1[j, 0], e1[j, 1], e1[j, 2], e2[j, 0], e2[j, 1], e2[j, 2],
                        nrm[j, 0], nrm[j, 1], nrm[j, 2], rad[j], aux[j])


@njit(cache=True, nogil=True)
def simulate(kind, pos, nrm, e1, e2, rad, aux, rb, L, ncell, stop_when_all_dirs):
    """Insert particles in order. Return (N_any, N_x, N_y, N_z), -1 if missing.

    kind: 0 disc, 1 cylinder, 2 square, 3 triangle, 4 spheroid.
    stop_when_all_dirs: 1 keeps inserting until x, y and z have each wrapped
    or the arrays run out; 0 stops at the first wrapping event.
    """
    n_insert = pos.shape[0]
    parent = np.empty(n_insert, np.int64)
    rank = np.empty(n_insert, np.int64)
    shift = np.zeros((n_insert, 3), np.int64)
    ncell3 = ncell * ncell * ncell
    head = np.empty(ncell3, np.int64)
    for c in range(ncell3):
        head[c] = -1
    nxt = np.empty(n_insert, np.int64)
    cell_edge = L / ncell
    Nx = -1
    Ny = -1
    Nz = -1
    Nany = -1
    rb_run = 0.0
    for i in range(n_insert):
        parent[i] = i
        rank[i] = 0
        shift[i, 0] = 0
        shift[i, 1] = 0
        shift[i, 2] = 0
        if rb[i] > rb_run:
            rb_run = rb[i]
        cx = _cell_of(pos[i, 0], L, ncell)
        cy = _cell_of(pos[i, 1], L, ncell)
        cz = _cell_of(pos[i, 2], L, ncell)
        cid = (cx * ncell + cy) * ncell + cz
        nxt[i] = head[cid]
        head[cid] = i
        r_int = rb[i] + rb_run + 1e-8
        if cell_edge > 0.0:
            rst = int(math.floor(r_int / cell_edge + 1e-12)) + 1
        else:
            rst = ncell
        # Visit each neighbouring cell once.
        if 2 * rst + 1 >= ncell:
            for nz in range(ncell):
                for ny in range(ncell):
                    for nx in range(ncell):
                        j = head[(nx * ncell + ny) * ncell + nz]
                        while j != -1:
                            if j < i:
                                s0, s1, s2 = _min_image_shift(
                                    pos[i, 0], pos[i, 1], pos[i, 2],
                                    pos[j, 0], pos[j, 1], pos[j, 2], L)
                                qjx = pos[j, 0] + s0 * L
                                qjy = pos[j, 1] + s1 * L
                                qjz = pos[j, 2] + s2 * L
                                ddx = qjx - pos[i, 0]
                                ddy = qjy - pos[i, 1]
                                ddz = qjz - pos[i, 2]
                                lim = rb[i] + rb[j] + 1e-10
                                if ddx * ddx + ddy * ddy + ddz * ddz <= lim * lim:
                                    if _contact(kind, i, j, qjx, qjy, qjz, pos, nrm, e1, e2, rad, aux):
                                        w0, w1, w2 = _union(i, j, s0, s1, s2, parent, shift, rank)
                                        if w0 != 0 and Nx < 0:
                                            Nx = i + 1
                                        if w1 != 0 and Ny < 0:
                                            Ny = i + 1
                                        if w2 != 0 and Nz < 0:
                                            Nz = i + 1
                                        if Nany < 0 and (w0 != 0 or w1 != 0 or w2 != 0):
                                            Nany = i + 1
                            j = nxt[j]
        else:
            for dz in range(-rst, rst + 1):
                nz = (cz + dz) % ncell
                for dy in range(-rst, rst + 1):
                    ny = (cy + dy) % ncell
                    for dx in range(-rst, rst + 1):
                        nx = (cx + dx) % ncell
                        j = head[(nx * ncell + ny) * ncell + nz]
                        while j != -1:
                            if j < i:
                                s0, s1, s2 = _min_image_shift(
                                    pos[i, 0], pos[i, 1], pos[i, 2],
                                    pos[j, 0], pos[j, 1], pos[j, 2], L)
                                qjx = pos[j, 0] + s0 * L
                                qjy = pos[j, 1] + s1 * L
                                qjz = pos[j, 2] + s2 * L
                                ddx = qjx - pos[i, 0]
                                ddy = qjy - pos[i, 1]
                                ddz = qjz - pos[i, 2]
                                lim = rb[i] + rb[j] + 1e-10
                                if ddx * ddx + ddy * ddy + ddz * ddz <= lim * lim:
                                    if _contact(kind, i, j, qjx, qjy, qjz, pos, nrm, e1, e2, rad, aux):
                                        w0, w1, w2 = _union(i, j, s0, s1, s2, parent, shift, rank)
                                        if w0 != 0 and Nx < 0:
                                            Nx = i + 1
                                        if w1 != 0 and Ny < 0:
                                            Ny = i + 1
                                        if w2 != 0 and Nz < 0:
                                            Nz = i + 1
                                        if Nany < 0 and (w0 != 0 or w1 != 0 or w2 != 0):
                                            Nany = i + 1
                            j = nxt[j]
        if Nany >= 0:
            if stop_when_all_dirs == 0:
                break
            if Nx >= 0 and Ny >= 0 and Nz >= 0:
                break
    out = np.empty(4, np.int64)
    out[0] = Nany
    out[1] = Nx
    out[2] = Ny
    out[3] = Nz
    return out


def choose_ncell(L, r_int_typical):
    """Grid with more than one cell when the box allows it.

    Correctness does not depend on this choice: the stencil radius is derived
    from the interaction range. The target edge is a bit under two radii so
    unit-particle runs keep a 3x3x3 stencil.
    """
    target = max(float(r_int_typical), L / 48.0)
    ncell = int(math.floor(L / target))
    if ncell < 1:
        ncell = 1
    if ncell > 48:
        ncell = 48
    return ncell


def frame_from_normal(n, spin):
    """Right-handed frame with the given normal and an in-plane spin angle."""
    n = np.asarray(n, dtype=np.float64)
    n = n / np.linalg.norm(n)
    if abs(n[2]) < 0.9:
        ref = np.array([0.0, 0.0, 1.0])
    else:
        ref = np.array([1.0, 0.0, 0.0])
    e1 = np.cross(ref, n)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    c = math.cos(spin)
    s = math.sin(spin)
    e1r = c * e1 + s * e2
    e2r = np.cross(n, e1r)
    return e1r, e2r, n


def random_normal(rng):
    z = rng.uniform(-1.0, 1.0)
    phi = rng.uniform(0.0, 2.0 * math.pi)
    r = math.sqrt(max(0.0, 1.0 - z * z))
    return np.array([r * math.cos(phi), r * math.sin(phi), z], dtype=np.float64)


def fisher_cosine(u, kappa):
    """Inverse CDF given in the specification. u in [0, 1]."""
    ek = math.exp(-2.0 * kappa)
    # u == 0 makes the log argument ek, which is fine in float64 for this kappa.
    arg = u + (1.0 - u) * ek
    if arg < 1e-300:
        arg = 1e-300
    c = 1.0 + math.log(arg) / kappa
    if c > 1.0:
        c = 1.0
    if c < -1.0:
        c = -1.0
    return c


def fisher_normal(rng, director, kappa):
    c = fisher_cosine(rng.random(), kappa)
    phi = rng.uniform(0.0, 2.0 * math.pi)
    s = math.sqrt(max(0.0, 1.0 - c * c))
    d = np.asarray(director, dtype=np.float64)
    d = d / np.linalg.norm(d)
    e1, e2, _ = frame_from_normal(d, 0.0)
    return s * math.cos(phi) * e1 + s * math.sin(phi) * e2 + c * d
