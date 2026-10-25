"""Finite-size fit of the single-axis winding estimator for thin discs (data: axes_L{16,24,32,64,96}.json from axes_check.py).
Model: eta_axis(L) = eta_inf + a * L**(-1/nu). Fit (i) with nu fixed at the 3D percolation value 0.8765 and (ii) with the exponent free.
The any-direction estimator is also listed. Run from workdir: python3 verification/axes_fss_fit.py"""
import json, math, numpy as np
from scipy.optimize import curve_fit
K = 4*math.pi/3; Ls = [16, 24, 32, 64, 96]; ax, axe, anyv, anye = [], [], [], []
for L in Ls:
    a = np.array(json.load(open(f"verification/axes_L{L}.json")), float)*K/L**3
    m = a[:, 1:].mean(1); ax.append(m.mean()); axe.append(m.std(ddof=1)/math.sqrt(len(m))); anyv.append(a[:, 0].mean()); anye.append(a[:, 0].std(ddof=1)/math.sqrt(len(a)))
Ls = np.array(Ls, float); ax = np.array(ax); axe = np.array(axe)
print("L, mean-of-axes eta (se), any-direction eta (se)")
for i, L in enumerate(Ls): print(f"{L:4.0f}  {ax[i]:.4f} ({axe[i]:.4f})   {anyv[i]:.4f} ({anye[i]:.4f})")
nu = 0.8765
f1 = lambda L, e, a: e + a*L**(-1/nu)
p, c = curve_fit(f1, Ls, ax, p0=[0.97, 1.0], sigma=axe, absolute_sigma=True); chi = (((ax-f1(Ls, *p))/axe)**2).sum()
print(f"fixed 1/nu = {1/nu:.3f}: eta_inf = {p[0]:.4f} +- {math.sqrt(c[0,0]):.4f}, a = {p[1]:.3f}, chi2 = {chi:.2f} for {len(Ls)-2} dof")
f2 = lambda L, e, a, b: e + a*L**(-b)
p2, c2 = curve_fit(f2, Ls, ax, p0=[0.97, 1.0, 1.1], sigma=axe, absolute_sigma=True, maxfev=20000); chi2 = (((ax-f2(Ls, *p2))/axe)**2).sum()
print(f"free exponent: eta_inf = {p2[0]:.4f} +- {math.sqrt(c2[0,0]):.4f}, a = {p2[1]:.3f}, b = {p2[2]:.2f} +- {math.sqrt(c2[2,2]):.2f}, chi2 = {chi2:.2f} for {len(Ls)-3} dof")
w = 1/np.array(anye)**2; big = Ls >= 24; print(f"any-direction plateau (L>=24, inverse-variance): {(w[big]*np.array(anyv)[big]).sum()/w[big].sum():.4f} +- {1/math.sqrt(w[big].sum()):.4f}")
print(f"published thin-disc value 0.9614; single-axis limit / 0.9614 = {p[0]/0.9614:.4f}")
