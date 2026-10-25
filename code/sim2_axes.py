"""Variant of sim2.first_perc (modes 0 and 1 only) that keeps inserting after the first wrap and records, per realisation,
N_any (first non-zero winding), and N_x, N_y, N_z = first insertion count at which a closing path has a non-zero winding component along that axis.
Used to test whether the 'any direction' criterion explains part of the offset from the published thin-disc threshold."""
import math, numpy as np
from numba import njit
from sim2 import disc_contact, gjk_cyl, find, random_rotations

@njit(cache=True)
def first_perc_axes(mode, pos, rot, rad, hth, L, ncell, brad):
    N = pos.shape[0]; cs = L/ncell
    head = -np.ones((ncell, ncell, ncell), np.int64); nxt = -np.ones(N, np.int64)
    par = np.arange(N); off = np.zeros((N, 3), np.int64); size = np.ones(N, np.int64); zero = np.zeros(3)
    ni = np.empty(3); nj = np.empty(3); res = -np.ones(4, np.int64)   # any, x, y, z
    for i in range(N):
        cx = int(pos[i, 0]/cs) % ncell; cy = int(pos[i, 1]/cs) % ncell; cz = int(pos[i, 2]/cs) % ncell
        for q in range(3): ni[q] = rot[i, q, 2]
        for dx in range(-1, 2):
            for dy in range(-1, 2):
                for dz in range(-1, 2):
                    j = head[(cx+dx) % ncell, (cy+dy) % ncell, (cz+dz) % ncell]
                    while j != -1:
                        dv0 = pos[j]-pos[i]; sh = -np.round(dv0/L); dv = dv0+sh*L; rr = brad[i]+brad[j]
                        if dv[0]**2+dv[1]**2+dv[2]**2 < rr*rr:
                            for q in range(3): nj[q] = rot[j, q, 2]
                            if mode == 0: hit = disc_contact(ni, rad[i], dv, nj, rad[j])
                            else: hit = gjk_cyl(zero, ni, rad[i], hth[i], dv, nj, rad[j], hth[j])
                            if hit:
                                shi = sh.astype(np.int64); ri, oi = find(par, off, i); rj, oj = find(par, off, j)
                                if ri == rj:
                                    w = oi+shi-oj
                                    if w[0] != 0 or w[1] != 0 or w[2] != 0:
                                        if res[0] < 0: res[0] = i+1
                                        for a in range(3):
                                            if w[a] != 0 and res[a+1] < 0: res[a+1] = i+1
                                        if res[1] >= 0 and res[2] >= 0 and res[3] >= 0: return res
                                else:
                                    if size[ri] >= size[rj]:
                                        par[rj] = ri; off[rj] = oi+shi-oj; size[ri] += size[rj]
                                    else:
                                        par[ri] = rj; off[ri] = oj-shi-oi; size[rj] += size[ri]
                        j = nxt[j]
        nxt[i] = head[cx, cy, cz]; head[cx, cy, cz] = i
    return res
