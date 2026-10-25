"""Extended continuum-percolation simulator for platelet fillers (periodic box, incremental insertion,
weighted union-find that detects the first wrapping cycle).

Contact modes
  0 : zero-thickness flat discs (exact analytic test)
  1 : finite-thickness cylinders (exact, GJK with analytic support functions)
  2 : sheets given as parametric surfaces (flat / cylindrical ripple / spherical cap), mesh-vs-implicit-surface test
Units: lengths in units of the mean disc radius.
"""
import math
import numpy as np
from numba import njit


# ----------------------------------------------------------------------------- small vector helpers
@njit(cache=True)
def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


@njit(cache=True)
def cross(a, b):
    out = np.empty(3)
    out[0] = a[1] * b[2] - a[2] * b[1]
    out[1] = a[2] * b[0] - a[0] * b[2]
    out[2] = a[0] * b[1] - a[1] * b[0]
    return out


# ----------------------------------------------------------------------------- thin-disc analytic contact
@njit(cache=True)
def disc_contact(ni, Ri, d, nj, Rj):
    """Exact test for two zero-thickness discs; disc i at origin, disc j at d."""
    m = cross(ni, nj)
    ms = math.sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2])
    if ms < 1e-9:
        h = abs(dot(d, ni))
        if h > 1e-12:
            return False
        dl = d - dot(d, ni) * ni
        return math.sqrt(dot(dl, dl)) < Ri + Rj
    u = m / ms
    c = dot(ni, nj)
    h2 = dot(d, nj)
    den = 1.0 - c * c
    p0 = (-h2 * c * ni + h2 * nj) / den
    b1 = dot(p0, u)
    q1 = dot(p0, p0)
    disc1 = b1 * b1 - (q1 - Ri * Ri)
    if disc1 < 0:
        return False
    q = p0 - d
    b2 = dot(q, u)
    q2 = dot(q, q)
    disc2 = b2 * b2 - (q2 - Rj * Rj)
    if disc2 < 0:
        return False
    s1 = math.sqrt(disc1)
    s2 = math.sqrt(disc2)
    return max(-b1 - s1, -b2 - s2) <= min(-b1 + s1, -b2 + s2)


# ----------------------------------------------------------------------------- GJK for cylinders
@njit(cache=True)
def sup_cyl(c, n, R, h, d, out):
    dn = dot(d, n)
    px = d[0] - dn * n[0]
    py = d[1] - dn * n[1]
    pz = d[2] - dn * n[2]
    l = math.sqrt(px * px + py * py + pz * pz)
    s = h if dn >= 0 else -h
    f = R / l if l > 1e-14 else 0.0
    out[0] = c[0] + s * n[0] + f * px
    out[1] = c[1] + s * n[1] + f * py
    out[2] = c[2] + s * n[2] + f * pz


@njit(cache=True)
def simplex_line(S, d):
    A = S[1].copy()
    B = S[0].copy()
    AB = B - A
    AO = -A
    if dot(AB, AO) > 0:
        d[:] = AO * dot(AB, AB) - AB * dot(AB, AO)
        return 2
    S[0] = A
    d[:] = AO
    return 1


@njit(cache=True)
def simplex_tri(S, d):
    A = S[2].copy()
    B = S[1].copy()
    C = S[0].copy()
    AB = B - A
    AC = C - A
    AO = -A
    ABC = cross(AB, AC)
    if dot(cross(ABC, AC), AO) > 0:
        if dot(AC, AO) > 0:
            S[0] = C
            S[1] = A
            d[:] = AO * dot(AC, AC) - AC * dot(AC, AO)
            return 2
        S[0] = B
        S[1] = A
        return simplex_line(S, d)
    if dot(cross(AB, ABC), AO) > 0:
        S[0] = B
        S[1] = A
        return simplex_line(S, d)
    if dot(ABC, AO) > 0:
        S[0] = C
        S[1] = B
        S[2] = A
        d[:] = ABC
        return 3
    S[0] = B
    S[1] = C
    S[2] = A
    d[:] = -ABC
    return 3


@njit(cache=True)
def simplex_tet(S, d):
    A = S[3].copy()
    B = S[2].copy()
    C = S[1].copy()
    D = S[0].copy()
    AB = B - A
    AC = C - A
    AD = D - A
    AO = -A
    ABC = cross(AB, AC)
    ACD = cross(AC, AD)
    ADB = cross(AD, AB)
    if dot(ABC, AO) > 0:
        S[0] = C
        S[1] = B
        S[2] = A
        return simplex_tri(S, d)
    if dot(ACD, AO) > 0:
        S[0] = D
        S[1] = C
        S[2] = A
        return simplex_tri(S, d)
    if dot(ADB, AO) > 0:
        S[0] = B
        S[1] = D
        S[2] = A
        return simplex_tri(S, d)
    return -1


@njit(cache=True)
def gjk_cyl(cA, nA, RA, hA, cB, nB, RB, hB):
    """True if two cylinders (centre, axis, radius, half-thickness) overlap."""
    S = np.zeros((4, 3))
    a = np.empty(3)
    b = np.empty(3)
    d = cB - cA
    if dot(d, d) < 1e-20:
        return True
    sup_cyl(cA, nA, RA, hA, d, a)
    sup_cyl(cB, nB, RB, hB, -d, b)
    S[0] = a - b
    d = -S[0]
    n = 1
    for it in range(60):
        if dot(d, d) < 1e-22:
            return True
        sup_cyl(cA, nA, RA, hA, d, a)
        sup_cyl(cB, nB, RB, hB, -d, b)
        p = a - b
        if dot(p, d) < 0:
            return False
        S[n] = p
        n += 1
        if n == 2:
            n = simplex_line(S, d)
        elif n == 3:
            n = simplex_tri(S, d)
        else:
            n = simplex_tet(S, d)
            if n == -1:
                return True
    return False


# ----------------------------------------------------------------------------- parametric-surface sheets
@njit(cache=True)
def surf_g(p, stype, sp0, sp1, sp2):
    """Implicit function whose zero set is the sheet surface (local frame)."""
    if stype == 0:
        return p[2]
    if stype == 1:  # ripple z = A sin(k x);   sp0 = A, sp1 = k
        return p[2] - sp0 * math.sin(sp1 * p[0])
    # spherical cap: sphere of radius Rc centred at (0,0,Rc); sp0 = Rc
    dx = p[0]
    dy = p[1]
    dz = p[2] - sp0
    return math.sqrt(dx * dx + dy * dy + dz * dz) - sp0


@njit(cache=True)
def in_domain(p, stype, sp0, sp1, sp2, R0, xtab, utab):
    """Is the surface point p inside the sheet's material disc?"""
    if stype == 0:
        return p[0] * p[0] + p[1] * p[1] <= R0 * R0
    if stype == 1:
        x = p[0]
        if x < xtab[0] or x > xtab[xtab.shape[0] - 1]:
            return False
        u = np.interp(x, xtab, utab)
        return u * u + p[1] * p[1] <= R0 * R0
    # cap: sp2 = maximal polar angle's z limit
    return p[2] <= sp2


@njit(cache=True)
def sheet_cross(Y, E, stype, sp0, sp1, sp2, R0, xtab, utab, g):
    m = Y.shape[0]
    for k in range(m):
        g[k] = surf_g(Y[k], stype, sp0, sp1, sp2)
    p = np.empty(3)
    for e in range(E.shape[0]):
        a = E[e, 0]
        b = E[e, 1]
        if g[a] * g[b] < 0:
            t = g[a] / (g[a] - g[b])
            for q in range(3):
                p[q] = Y[a, q] + t * (Y[b, q] - Y[a, q])
            if in_domain(p, stype, sp0, sp1, sp2, R0, xtab, utab):
                return True
    return False


@njit(cache=True)
def sheet_contact(RiM, dv, RjM, V, E, stype, sp0, sp1, sp2, R0, xtab, utab, Y, g):
    """Two parametric sheets: i at origin (rotation RiM), j at dv (rotation RjM)."""
    M = np.zeros((3, 3))
    Mt = np.zeros((3, 3))
    b = np.zeros(3)
    b2 = np.zeros(3)
    for q in range(3):
        for r in range(3):
            for k in range(3):
                M[q, r] += RjM[k, q] * RiM[k, r]
                Mt[q, r] += RiM[k, q] * RjM[k, r]
        for k in range(3):
            b[q] -= RjM[k, q] * dv[k]
            b2[q] += RiM[k, q] * dv[k]
    m = V.shape[0]
    for k in range(m):
        for q in range(3):
            Y[k, q] = M[q, 0] * V[k, 0] + M[q, 1] * V[k, 1] + M[q, 2] * V[k, 2] + b[q]
    if sheet_cross(Y[:m], E, stype, sp0, sp1, sp2, R0, xtab, utab, g):
        return True
    for k in range(m):
        for q in range(3):
            Y[k, q] = Mt[q, 0] * V[k, 0] + Mt[q, 1] * V[k, 1] + Mt[q, 2] * V[k, 2] + b2[q]
    return sheet_cross(Y[:m], E, stype, sp0, sp1, sp2, R0, xtab, utab, g)


# ----------------------------------------------------------------------------- union-find with periodic offsets
@njit(cache=True)
def find(par, off, x):
    o = np.zeros(3, np.int64)
    while par[x] != x:
        o += off[x]
        x = par[x]
    return x, o


@njit(cache=True)
def first_perc(mode, pos, rot, rad, hth, L, ncell, brad, V, E, stype, sp0, sp1, sp2, R0, xtab, utab):
    """Insert objects 0..N-1 in order; return the number inserted when the first wrapping cluster appears (-1 if never)."""
    N = pos.shape[0]
    cs = L / ncell
    head = -np.ones((ncell, ncell, ncell), np.int64)
    nxt = -np.ones(N, np.int64)
    par = np.arange(N)
    off = np.zeros((N, 3), np.int64)
    size = np.ones(N, np.int64)
    zero = np.zeros(3)
    Y = np.empty((V.shape[0], 3))
    g = np.empty(V.shape[0])
    ni = np.empty(3)
    nj = np.empty(3)
    for i in range(N):
        cx = int(pos[i, 0] / cs) % ncell
        cy = int(pos[i, 1] / cs) % ncell
        cz = int(pos[i, 2] / cs) % ncell
        for q in range(3):
            ni[q] = rot[i, q, 2]
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                for dz in range(-1, 2):
                    j = head[(cx + dx) % ncell, (cy + dy) % ncell, (cz + dz) % ncell]
                    while j != -1:
                        dv0 = pos[j] - pos[i]
                        sh = -np.round(dv0 / L)
                        dv = dv0 + sh * L
                        rr = brad[i] + brad[j]
                        if dv[0] ** 2 + dv[1] ** 2 + dv[2] ** 2 < rr * rr:
                            for q in range(3):
                                nj[q] = rot[j, q, 2]
                            hit = False
                            if mode == 0:
                                hit = disc_contact(ni, rad[i], dv, nj, rad[j])
                            elif mode == 1:
                                hit = gjk_cyl(zero, ni, rad[i], hth[i], dv, nj, rad[j], hth[j])
                            else:
                                hit = sheet_contact(rot[i], dv, rot[j], V, E, stype, sp0, sp1, sp2, R0, xtab, utab, Y, g)
                            if hit:
                                shi = sh.astype(np.int64)
                                ri, oi = find(par, off, i)
                                rj, oj = find(par, off, j)
                                if ri == rj:
                                    w = oi + shi - oj
                                    if w[0] != 0 or w[1] != 0 or w[2] != 0:
                                        return i + 1
                                else:
                                    if size[ri] >= size[rj]:
                                        par[rj] = ri
                                        off[rj] = oi + shi - oj
                                        size[ri] += size[rj]
                                    else:
                                        par[ri] = rj
                                        off[ri] = oj - shi - oi
                                        size[rj] += size[ri]
                        j = nxt[j]
        nxt[i] = head[cx, cy, cz]
        head[cx, cy, cz] = i
    return -1


# ----------------------------------------------------------------------------- python-side helpers
def random_rotations(n, rng):
    q = rng.normal(size=(n, 4))
    q /= np.linalg.norm(q, axis=1)[:, None]
    w, x, y, z = q.T
    R = np.empty((n, 3, 3))
    R[:, 0, 0] = 1 - 2 * (y * y + z * z)
    R[:, 0, 1] = 2 * (x * y - z * w)
    R[:, 0, 2] = 2 * (x * z + y * w)
    R[:, 1, 0] = 2 * (x * y + z * w)
    R[:, 1, 1] = 1 - 2 * (x * x + z * z)
    R[:, 1, 2] = 2 * (y * z - x * w)
    R[:, 2, 0] = 2 * (x * z - y * w)
    R[:, 2, 1] = 2 * (y * z + x * w)
    R[:, 2, 2] = 1 - 2 * (x * x + y * y)
    return R


def rotation_with_normal(nrm, rng):
    """Rotation matrices whose third column is the given unit normal, random in-plane angle."""
    n = nrm.shape[0]
    a = np.where(np.abs(nrm[:, [0]]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    e1 = np.cross(nrm, a)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(nrm, e1)
    ang = rng.random(n) * 2 * np.pi
    c, s = np.cos(ang)[:, None], np.sin(ang)[:, None]
    f1 = c * e1 + s * e2
    f2 = -s * e1 + c * e2
    R = np.empty((n, 3, 3))
    R[:, :, 0] = f1
    R[:, :, 1] = f2
    R[:, :, 2] = nrm
    return R


EMPTY = dict(V=np.zeros((1, 3)), E=np.zeros((1, 2), np.int64), stype=0, sp0=0.0, sp1=0.0, sp2=0.0, R0=1.0,
             xtab=np.zeros(2), utab=np.zeros(2))


def make_sheet(kind, R0=1.0, h=None, amp=0.0, lam=1.0, Rc=1.0, nrim=48):
    """Mesh + implicit-surface description of a sheet whose material disc has radius R0.
    kind: 'flat', 'ripple' (isometric cylindrical corrugation z=A sin kx), 'cap' (isometric spherical cap)."""
    if h is None:
        h = R0 / 8.0
    ng = int(math.ceil(R0 / h))
    pts = []
    idx = {}
    for i in range(-ng, ng + 1):
        for j in range(-ng, ng + 1):
            u, v = i * h, j * h
            if u * u + v * v <= (R0 * 0.999) ** 2:
                idx[(i, j)] = len(pts)
                pts.append((u, v))
    edges = []
    for (i, j), a in idx.items():
        for (di, dj) in ((1, 0), (0, 1)):
            b = idx.get((i + di, j + dj))
            if b is not None:
                edges.append((a, b))
    nint = len(pts)
    for k in range(nrim):
        t = 2 * math.pi * k / nrim
        pts.append((R0 * math.cos(t), R0 * math.sin(t)))
    for k in range(nrim):
        edges.append((nint + k, nint + (k + 1) % nrim))
    P = np.array(pts)
    d = dict(E=np.array(edges, np.int64), R0=R0)
    if kind == 'flat':
        d.update(V=np.column_stack([P[:, 0], P[:, 1], np.zeros(len(P))]), stype=0, sp0=0.0, sp1=0.0, sp2=0.0,
                 xtab=np.zeros(2), utab=np.zeros(2))
    elif kind == 'ripple':
        k = 2 * math.pi / lam
        xs = np.linspace(-1.5 * R0, 1.5 * R0, 6001)
        us = np.concatenate([[0], np.cumsum(0.5 * (np.sqrt(1 + (amp * k * np.cos(k * xs[:-1])) ** 2) +
                                              np.sqrt(1 + (amp * k * np.cos(k * xs[1:])) ** 2)) * np.diff(xs))])
        us = us - np.interp(0.0, xs, us)  # u(0)=0
        x_of_u = np.interp(P[:, 0], us, xs)
        d.update(V=np.column_stack([x_of_u, P[:, 1], amp * np.sin(k * x_of_u)]), stype=1, sp0=amp, sp1=k, sp2=0.0,
                 xtab=xs, utab=us)
    elif kind == 'cap':
        rho = np.hypot(P[:, 0], P[:, 1])
        psi = np.arctan2(P[:, 1], P[:, 0])
        beta = rho / Rc
        d.update(V=np.column_stack([Rc * np.sin(beta) * np.cos(psi), Rc * np.sin(beta) * np.sin(psi),
                                    Rc * (1 - np.cos(beta))]), stype=2, sp0=Rc, sp1=0.0,
                 sp2=Rc * (1 - math.cos(R0 / Rc)), xtab=np.zeros(2), utab=np.zeros(2))
    else:
        raise ValueError(kind)
    d['brad'] = float(np.linalg.norm(d['V'], axis=1).max()) + 1e-9
    return d


def run_once(mode, pos, rot, rad, hth, L, brad, sheet=None, ncell=None):
    if L < 4.0 * brad.max():
        raise ValueError("box side L must be at least 4 x the largest bounding radius: contacts are recorded on the minimum image only")
    s = sheet if sheet is not None else EMPTY
    if ncell is None:
        ncell = max(3, int(L // (2 * brad.max())))
    return first_perc(mode, pos, rot, rad, hth, L, ncell, brad, s['V'], s['E'], s['stype'], s['sp0'], s['sp1'],
                      s['sp2'], s['R0'], s['xtab'], s['utab'])
