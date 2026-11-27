"""Figure support: configurations from the production simulator (sim4.py), with the spanning cluster and growth history.
Run from this folder; sim4.py, align.py and sim2-style helpers are imported from ../../code (see fig_common.py)."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from numba import njit
import sim4
from sim4 import (disc_contact, gjk_cyl, sheet_contact, find, run_once, random_rotations, rotation_with_normal,
                  make_sheet, EMPTY)

@njit(cache=True)
def trace_perc(mode, stop, pos, rot, rad, hth, L, ncell, brad, V, E, stype, sp0, sp1, sp2, R0, xtab, utab, offs):
    """Same insertion and contact logic as sim4.first_perc.  Inserts objects in order until the first wrapping cycle (or until
    `stop` objects if stop > 0).  Returns (count, root label of each inserted object, size of the largest cluster after each insertion,
    root of the wrapping cluster or -1)."""
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
    hist = np.zeros(N, np.int64)
    big = 1
    lim = N
    if stop > 0:
        lim = stop
    count = -1
    wroot = -1
    for i in range(lim):
        cx = int(pos[i, 0] / cs) % ncell
        cy = int(pos[i, 1] / cs) % ncell
        cz = int(pos[i, 2] / cs) % ncell
        for q in range(3):
            ni[q] = rot[i, q, 2]
        done = False
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
                                hit = sheet_contact(rot[i], dv, rot[j], V, E, stype, sp0, sp1, sp2, R0, xtab, utab, Y, g, offs, rad[i], rad[j])
                            if hit:
                                shi = sh.astype(np.int64)
                                ri, oi = find(par, off, i)
                                rj, oj = find(par, off, j)
                                if ri == rj:
                                    w = oi + shi - oj
                                    if w[0] != 0 or w[1] != 0 or w[2] != 0:
                                        count = i + 1
                                        wroot = ri
                                        done = True
                                else:
                                    if size[ri] >= size[rj]:
                                        par[rj] = ri
                                        off[rj] = oi + shi - oj
                                        size[ri] += size[rj]
                                        if size[ri] > big:
                                            big = size[ri]
                                    else:
                                        par[ri] = rj
                                        off[ri] = oj - shi - oi
                                        size[rj] += size[ri]
                                        if size[rj] > big:
                                            big = size[rj]
                        if done:
                            break
                        j = nxt[j]
                    if done:
                        break
                if done:
                    break
            if done:
                break
        hist[i] = big
        if done:
            n_in = i + 1
            lab = np.empty(n_in, np.int64)
            for k in range(n_in):
                r, o = find(par, off, k)
                lab[k] = r
            return count, lab, hist[:n_in], wroot
        nxt[i] = head[cx, cy, cz]
        head[cx, cy, cz] = i
    n_in = lim
    lab = np.empty(n_in, np.int64)
    for k in range(n_in):
        r, o = find(par, off, k)
        lab[k] = r
    return -1, lab, hist[:n_in], -1

def run_trace(cfg, stop=0):
    s = cfg.get('sheet') or EMPTY
    L = cfg['L']; brad = cfg['brad']
    ncell = max(3, int(L // (2 * brad.max())))
    return trace_perc(cfg['mode'], stop, cfg['pos'], cfg['rot'], cfg['rad'], cfg['hth'], L, ncell, brad, s['V'], s['E'], s['stype'],
                      s['sp0'], s['sp1'], s['sp2'], s['R0'], s['xtab'], s['utab'], np.zeros(1))

def pair_order(pos, nrm, L, rmax=2.0):
    from align import shell_P2
    return shell_P2(pos, nrm, L, rmax)[0]

def make_cfg(kind, L=16.0, seed=1, A=30.0, N=1, sig=0.0, ell=4.0, kappa=0.0, amp_slope=0.0, lam=1.0, nfac=2.5, rcap=6.0):
    """kind: 'thin' (zero thickness discs), 'cyl' (finite-thickness cylinders, with optional tactoid layer number N, log-normal sigma,
    domain alignment ell/kappa) or 'ripple' (rippled zero-thickness sheets of peak slope amp_slope and wavelength lam)."""
    from align import aligned_normals
    rng = np.random.default_rng(seed)
    n = int(nfac * L ** 3)
    if kind == 'ripple':
        amp = amp_slope * lam / (2 * math.pi)
        sh = make_sheet('ripple', 1.0, h=1 / 8, amp=amp, lam=lam)
        pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng)
        rad = np.ones(n); brad = np.full(n, sh['brad'])
        return dict(kind=kind, mode=2, L=L, pos=pos, rot=rot, rad=rad, hth=np.zeros(n), brad=brad, sheet=sh, N=1, A=1e9,
                    amp=amp, lam=lam, kk=2 * math.pi / lam, slope=amp_slope)
    R = rng.lognormal(-0.5 * sig ** 2, sig, size=n) if sig > 0 else np.ones(n)
    R = R[R <= rcap]; n = len(R)
    pos = rng.random((n, 3)) * L
    if kappa > 0:
        nrm = aligned_normals(pos, L, ell, kappa, rng); rot = rotation_with_normal(nrm, rng)
    else:
        rot = random_rotations(n, rng)
    if kind == 'thin':
        hth = np.zeros(n); mode = 0; brad = R.copy()
    else:
        hth = np.full(n, N / A); mode = 1; brad = np.sqrt(R ** 2 + hth ** 2)
    return dict(kind=kind, mode=mode, L=L, pos=pos, rot=rot, rad=R, hth=hth, brad=brad, sheet=None, N=N, A=A, sig=sig)
