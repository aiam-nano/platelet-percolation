"""Checks of the zero-thickness polygon contact (mode 4): mean excluded volume of two randomly oriented plates
against the analytic V_ex = A*P/2 (A area, P perimeter; Isihara/Kihara for flat convex bodies), and a many-sided
polygon against the analytic disc test."""
import numpy as np, math
from numba import njit
from sim2e import *

def polygon(kind):
    """vertices (nv,3) of a unit-equivalent-radius plate (area = pi), centred at the centroid"""
    if kind == 'square':
        s = math.sqrt(math.pi); h = s / 2
        P = np.array([[-h, -h], [h, -h], [h, h], [-h, h]])
    elif kind == 'triangle':
        s = math.sqrt(4 * math.pi / math.sqrt(3)); R = s / math.sqrt(3)
        P = np.array([[R * math.cos(a), R * math.sin(a)] for a in (math.pi / 2, math.pi / 2 + 2 * math.pi / 3, math.pi / 2 + 4 * math.pi / 3)])
    else:   # n-gon of equal area
        n = int(kind); R = math.sqrt(2 * math.pi / (n * math.sin(2 * math.pi / n)))
        P = np.array([[R * math.cos(2 * math.pi * k / n), R * math.sin(2 * math.pi * k / n)] for k in range(n)])
    return np.column_stack([P, np.zeros(len(P))])

def perimeter(V):
    P = V[:, :2]; return float(np.sum(np.linalg.norm(np.roll(P, -1, axis=0) - P, axis=1)))

@njit(cache=True)
def pair_loop(Rots, dvs, V, rad):
    hit = 0; hitd = 0
    z = np.zeros(3)
    for k in range(Rots.shape[0] // 2):
        Ri = Rots[2 * k]; Rj = Rots[2 * k + 1]
        if poly_contact(Ri, dvs[k], Rj, V): hit += 1
        ni = np.empty(3); nj = np.empty(3)
        for q in range(3): ni[q] = Ri[q, 2]; nj[q] = Rj[q, 2]
        if disc_contact(ni, rad, dvs[k], nj, rad): hitd += 1
    return hit, hitd

if __name__ == "__main__":
    rng = np.random.default_rng(11)
    npair = 20_000_000
    for kind in ('square', 'triangle', '720'):
        V = polygon(kind); Rc = float(np.linalg.norm(V, axis=1).max()); box = 2 * Rc
        tot = 0; hd = 0; nb = 0
        for b in range(10):
            m = npair // 10
            Rots = random_rotations(2 * m, rng)
            dvs = (rng.random((m, 3)) - 0.5) * 2 * box      # cube of half-side box contains the contact region
            h, hdd = pair_loop(Rots, dvs, V, 1.0)
            tot += h; hd += hdd; nb += m
        vol = (2 * box) ** 3
        Vex = tot / nb * vol; err = Vex / math.sqrt(max(tot, 1))
        A = 0.5 * abs(np.sum(V[:, 0] * np.roll(V[:, 1], -1) - np.roll(V[:, 0], -1) * V[:, 1]))
        an = A * perimeter(V) / 2
        print(f"{kind:9s} V_ex(MC)={Vex:.4f}+-{err:.4f}  analytic A*P/2={an:.4f}   (same pairs, analytic disc of r=1: {hd / nb * vol:.4f}, pi^2={math.pi**2:.4f})", flush=True)
