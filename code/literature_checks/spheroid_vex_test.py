"""Mean excluded volume of two identical randomly oriented oblate spheroids against V_ex = 2V + S M / (2 pi) (Kihara), with
S the surface area and M the integral of the mean curvature of one spheroid (computed here by quadrature).
Tests gjk_ell (mode 3 contact) independently of any percolation run."""
import math, numpy as np
from numba import njit
from sim2e import gjk_ell, random_rotations
def SM(a, c, n=200001):
    th = np.linspace(1e-9, math.pi - 1e-9, n); rho = a * np.sin(th); z = c * np.cos(th)
    rp = a * np.cos(th); zp = -c * np.sin(th); rpp = -a * np.sin(th); zpp = -c * np.cos(th)
    sp = np.sqrt(rp ** 2 + zp ** 2); km = np.abs((rp * zpp - rpp * zp) / sp ** 3); kp = np.abs(zp / (rho * sp))
    dS = 2 * math.pi * rho * sp; S = np.trapezoid(dS, th); M = np.trapezoid(0.5 * (km + kp) * dS, th); return S, M
@njit(cache=True)
def count(Rots, dvs, a, c):
    hit = 0; z = np.zeros(3); ni = np.empty(3); nj = np.empty(3)
    for k in range(dvs.shape[0]):
        for q in range(3): ni[q] = Rots[2 * k, q, 2]; nj[q] = Rots[2 * k + 1, q, 2]
        if gjk_ell(z, ni, a, c, dvs[k], nj, a, c): hit += 1
    return hit
rng = np.random.default_rng(5)
for c in (1.0, 1 / 3, 0.1, 0.01):
    a = 1.0; rm = 2 * max(a, c); tot = 0; nb = 0
    for b in range(5):
        m = 2_000_000; R = random_rotations(2 * m, rng); d = (rng.random((m, 3)) - 0.5) * 2 * rm
        tot += count(R, d, a, c); nb += m
    vol = (2 * rm) ** 3; Vex = tot / nb * vol; err = Vex / math.sqrt(tot)
    S, M = SM(a, c); V = 4 / 3 * math.pi * a * a * c; an = 2 * V + S * M / (2 * math.pi)
    print(f"c/a={c:7.4f}: V_ex MC {Vex:.4f}+-{err:.4f}   Kihara 2V+SM/(2pi) {an:.4f}   ratio {Vex / an:.4f}", flush=True)
