"""Candidate explanations of the 1.6 per cent offset of the thin-disc threshold (25 Oct 2026, action A1).
Uses the stored first-wrap thresholds (fss_thin.json, fss_L64/96/128.json) and the per-axis arrays (axes_L*.json); nothing is simulated
except the tangent-contact test at the end.  Run from this folder:  python3 offset_diagnostics.py
Reference value eta_c = 0.9614 (Yi and Tawerghi 2009)."""
import json, math, numpy as np
from scipy.stats import poisson
from scipy.optimize import brentq, curve_fit
REF = 0.9614
data = {float(k): np.array(v) for k, v in json.load(open("fss_thin.json")).items()}
for L, f in ((64.0, "fss_L64.json"), (96.0, "fss_L96.json"), (128.0, "fss_L128.json")):
    data[L] = np.array(json.load(open(f)))
# 25 Oct night: 60 further realisations at L = 128 (seed 5128) and 40 at L = 192 (seed 5192), run with fss_big.py
data[128.0] = np.concatenate([data[128.0], np.array(json.load(open("fss_L128_extra.json")))])
data[192.0] = np.array(json.load(open("fss_L192.json")))
Ls = sorted(data)
print("1. Estimator: mean, median, and Poisson-ensemble 50 per cent crossing of the same sequential realisations")
print("   (the Poisson ensemble has a fluctuating particle number; P(wrap at mean density) = mean_i P(Poisson(N) >= N_i))")
res = {}
for L in Ls:
    eta = data[L]; f = 4 * math.pi / 3 / L**3
    Nc = eta / f
    P = lambda Nbar: np.mean(poisson.sf(np.round(Nc) - 1, Nbar)) - 0.5
    Nb = brentq(P, Nc.min() * 0.8, Nc.max() * 1.2)
    se = eta.std(ddof=1) / math.sqrt(len(eta))
    res[L] = (eta.mean(), se, np.median(eta), Nb * f)
    print(f"   L={L:5.0f} n={len(eta):4d}  mean {eta.mean():.4f}+-{se:.4f}  median {np.median(eta):.4f}  Poisson-crossing {Nb*f:.4f}   (sd of single realisation {eta.std(ddof=1):.3f}, Poisson sd of N/N {1/math.sqrt(Nc.mean()):.4f})")
big = [L for L in Ls if L >= 32]
w = np.array([1 / res[L][1] ** 2 for L in big])
for name, k in (("mean", 0), ("median", 2), ("Poisson crossing", 3)):
    v = np.array([res[L][k] for L in big]); print(f"   weighted over L>=32, {name:17s}: {np.sum(w*v)/np.sum(w):.4f}  ({np.sum(w*v)/np.sum(w)/REF-1:+.2%} against {REF})")
print("   Result: switching to a Poisson number of particles, or to the median, moves the plateau by well under 0.5 per cent; it cannot account for 1.6 per cent.")

print("\n2. Crossing of the wrapping-probability curves of two box sizes (no extrapolation in L)")
def Pcurve(eta, x):
    return np.searchsorted(np.sort(eta), x, side="right") / len(eta)
rng = np.random.default_rng(3)
grid = np.linspace(0.90, 1.05, 601)
def crossing(a, b):
    d = Pcurve(a, grid) - Pcurve(b, grid)
    # a larger box has the steeper curve: it is below 0.5 below eta_c and above it above eta_c, so look for the sign change of (large - small) near P ~ 0.5
    s = np.where(np.diff(np.sign(d)) != 0)[0]
    s = [i for i in s if 0.2 < Pcurve(a, grid[i]) < 0.8]
    return grid[s[0]] if s else np.nan
for (A, B) in ((16, 32), (24, 48), (32, 64), (48, 96), (32, 48), (48, 64)):
    base = crossing(data[float(B)], data[float(A)])
    bs = [crossing(rng.choice(data[float(B)], len(data[float(B)])), rng.choice(data[float(A)], len(data[float(A)]))) for _ in range(400)]
    bs = np.array([x for x in bs if x == x])
    print(f"   L={A:3d} and {B:3d}: crossing {base:.4f}  bootstrap {np.mean(bs):.4f} +- {np.std(bs):.4f}  ({len(bs)}/400 resolved)")
print("   Crossings are noisy (the curves are nearly parallel at these sizes) but sit at 0.97 to 0.98, not near 0.9614.")

print("\n3. Can a slow finite-size trend hide 0.9614?  Fit eta(L) = eta_inf + a L^-p to the means, L >= 12")
Lx = np.array([L for L in Ls if L >= 12]); y = np.array([res[L][0] for L in Lx]); s = np.array([res[L][1] for L in Lx])
def chi2(par, f): return float(np.sum(((y - f(Lx, *par)) / s) ** 2))
fr = lambda L, e, a, p: e + a * L ** (-p)
popt, pc = curve_fit(fr, Lx, y, p0=[0.975, 1, 1], sigma=s, absolute_sigma=True, maxfev=20000)
print(f"   free eta_inf, a, p: eta_inf {popt[0]:.4f}+-{math.sqrt(pc[0,0]):.4f}, p {popt[2]:.2f}, chi2 {chi2(popt, fr):.1f} for {len(Lx)-3} dof")
f2 = lambda L, a, p: REF + a * L ** (-p)
p2, _ = curve_fit(f2, Lx, y, p0=[1, 0.3], sigma=s, absolute_sigma=True, maxfev=20000)
print(f"   eta_inf fixed at {REF}: best a {p2[0]:.3f}, p {p2[1]:.2f}, chi2 {float(np.sum(((y-f2(Lx,*p2))/s)**2)):.1f} for {len(Lx)-2} dof")
f3 = lambda L, a: 0.9770 + a * L ** (-1.141)
p3, _ = curve_fit(f3, Lx, y, p0=[1], sigma=s, absolute_sigma=True)
print(f"   eta_inf fixed at 0.9770 with the 3D exponent 1/nu = 1.141: a {p3[0]:.3f}, chi2 {float(np.sum(((y-f3(Lx,*p3))/s)**2)):.1f} for {len(Lx)-1} dof")

print("\n4. Single-axis wrapping (the periodic analogue of spanning in one direction)")
for L in (16, 24, 96):
    a = np.array(json.load(open(f"axes_L{L}.json")), float) / L**3 * 4 * math.pi / 3
    print(f"   L={L:3d}: any {a[:,0].mean():.4f}  x {a[:,1].mean():.4f}  y {a[:,2].mean():.4f}  z {a[:,3].mean():.4f}")
print("   Extrapolation of the single-axis values in 1/L (axes_fss_fit.log) gives 0.9756 +- 0.0012, 1.5 per cent above 0.9614.")

print("\n5. Tangent contacts: discs with radius scaled by 1 +- 1e-9 give identical first-wrap indices")
import sys; sys.path.insert(0, "..")
from sim2 import random_rotations, run_once
rng = np.random.default_rng(11); L = 12.0; diff = 0; n_test = 40
for _ in range(n_test):
    n = int(0.35 * L**3); pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng)
    a = run_once(0, pos, rot, np.ones(n) * (1 - 1e-9), np.zeros(n), L, np.ones(n) * 1.001)
    b = run_once(0, pos, rot, np.ones(n) * (1 + 1e-9), np.zeros(n), L, np.ones(n) * 1.001)
    diff += (a != b)
print(f"   {diff} of {n_test} configurations differ.")
