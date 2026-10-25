"""Independent check of disc_contact: Monte Carlo mean excluded volume of two randomly oriented unit discs vs analytic pi^2."""
import sys, math, numpy as np
sys.path.insert(0, "..")
from numba import njit
from sim2 import disc_contact, random_rotations

@njit(cache=True)
def count(N, seed):
    np.random.seed(seed)
    hits = 0
    for _ in range(N):
        # random normals
        a = np.random.normal(0., 1., 3); a /= np.sqrt((a*a).sum())
        b = np.random.normal(0., 1., 3); b /= np.sqrt((b*b).sum())
        # centre uniform in ball radius 2
        while True:
            d = np.random.uniform(-2., 2., 3)
            if (d*d).sum() <= 4.0:
                break
        if disc_contact(a, 1.0, d, b, 1.0):
            hits += 1
    return hits

N = 40_000_000
h = count(N, 7)
p = h / N
Vball = 4/3*math.pi*8
vex = p*Vball; err = Vball*math.sqrt(p*(1-p)/N)
print(f"Vex_MC = {vex:.5f} +- {err:.5f}   analytic pi^2 = {math.pi**2:.5f}   ratio = {vex/math.pi**2:.5f}")
