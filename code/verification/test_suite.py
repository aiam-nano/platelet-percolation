"""Verification suite for the percolation simulator (sim2.py) and conventions.
Each test prints PASS/FAIL with the number behind the decision; results are written to test_report.json.
Run from this directory: python3 test_suite.py [quick]
Tolerances are stated in each test and were fixed before the run."""
import sys, json, math, time
import numpy as np
from numba import njit
sys.path.insert(0, "..")
from sim2 import disc_contact, gjk_cyl, random_rotations, rotation_with_normal, run_once
from align import aligned_normals, fisher
QUICK = len(sys.argv) > 1 and sys.argv[1] == "quick"
REP = []
def record(name, ok, detail):
    REP.append(dict(test=name, passed=bool(ok), detail=detail)); print(("PASS " if ok else "FAIL ") + name + " :: " + detail, flush=True)

rng = np.random.default_rng(2026)
def runit(k):
    v = rng.normal(size=(k, 3)); return v / np.linalg.norm(v, axis=1)[:, None]

# ---- T1: analytic disc contact against an independent brute-force piercing test (tolerance: <= 0.1 % mismatches)
def pierce(c1, n1, R1, c2, n2, R2, K=1440):
    a = np.cross(n1, [1, 0, 0])
    if np.linalg.norm(a) < 1e-3: a = np.cross(n1, [0, 1, 0])
    a /= np.linalg.norm(a); b = np.cross(n1, a)
    t = np.linspace(0, 2 * np.pi, K, endpoint=False)
    P = c1 + R1 * (np.outer(np.cos(t), a) + np.outer(np.sin(t), b)); Q = np.roll(P, -1, axis=0)
    s1 = (P - c2) @ n2; s2 = (Q - c2) @ n2; cr = (s1 * s2 < 0)
    if not cr.any(): return False
    w = s1[cr] / (s1[cr] - s2[cr]); X = P[cr] + w[:, None] * (Q[cr] - P[cr])
    return bool((np.linalg.norm(X - c2, axis=1) < R2).any())
n = 3000 if not QUICK else 800; bad = 0; pos = 0
for _ in range(n):
    c2 = rng.uniform(-2, 2, 3); n1 = runit(1)[0]; n2 = runit(1)[0]; r2 = rng.uniform(0.5, 1.5)
    A = disc_contact(n1, 1.0, c2, n2, r2); B = pierce(np.zeros(3), n1, 1.0, c2, n2, r2) or pierce(c2, n2, r2, np.zeros(3), n1, 1.0)
    bad += (A != B); pos += A
record("T1 disc_contact vs brute force", bad <= 0.001 * n, f"{bad} mismatches in {n} pairs ({pos} contacts)")

# ---- T2: GJK cylinder test in the thin limit against the analytic disc test (tolerance: 0 mismatches except sub-1e-6 tangencies)
bad = 0; n = 4000 if not QUICK else 1000; pos = 0
for _ in range(n):
    d = rng.uniform(-2.5, 2.5, 3); n1 = runit(1)[0]; n2 = runit(1)[0]; r2 = rng.uniform(0.6, 1.4)
    A = disc_contact(n1, 1.0, d, n2, r2); B = gjk_cyl(np.zeros(3), n1, 1.0, 1e-6, d, n2, r2, 1e-6)
    bad += (A != B); pos += A
record("T2 GJK thin limit vs analytic disc", bad <= 0.001 * n, f"{bad} mismatches in {n} pairs ({pos} contacts)")

# ---- T3: mean excluded volume (two-body Monte Carlo) against the Isihara formula V_ex = 2V + S M / (2 pi)
@njit(cache=True)
def mc_vex(Rr, h, box, N, seed):
    np.random.seed(seed); hits = 0
    c0 = np.zeros(3); n1 = np.empty(3); n2 = np.empty(3); d = np.empty(3)
    for _ in range(N):
        for q in range(3):
            d[q] = (np.random.random() * 2 - 1) * box
        v = np.random.normal(0, 1, 3); n1 = v / np.sqrt((v * v).sum())
        v = np.random.normal(0, 1, 3); n2 = v / np.sqrt((v * v).sum())
        if gjk_cyl(c0, n1, Rr, h, d, n2, Rr, h): hits += 1
    return hits
worst = 0.0; det = []
for t in (0.05, 0.4, 1.0, 2.0):
    Rr = 1.0; h = t / 2; V = math.pi * Rr * Rr * t; S = 2 * math.pi * Rr**2 + 2 * math.pi * Rr * t; M = math.pi * t + math.pi**2 * Rr
    Vex = 2 * V + S * M / (2 * math.pi); box = 2 * math.sqrt(Rr * Rr + h * h) + 0.01
    N = 4_000_000 if not QUICK else 1_000_000; hits = mc_vex(Rr, h, box, N, int(100 * t) + 7)
    p = hits / N; est = p * (2 * box)**3; se = (2 * box)**3 * math.sqrt(p * (1 - p) / N); z = (est - Vex) / se
    worst = max(worst, abs(z)); det.append(f"t={t}: MC {est:.4f}+-{se:.4f} vs {Vex:.4f} (z={z:+.1f})")
record("T3 excluded volume vs Isihara (cylinders)", worst < 4, "; ".join(det))
@njit(cache=True)
def mc_vex_disc(box, N, seed):
    np.random.seed(seed); hits = 0; d = np.empty(3); c0 = np.zeros(3)
    for _ in range(N):
        for q in range(3): d[q] = (np.random.random() * 2 - 1) * box
        v = np.random.normal(0, 1, 3); n1 = v / np.sqrt((v * v).sum()); v = np.random.normal(0, 1, 3); n2 = v / np.sqrt((v * v).sum())
        if disc_contact(n1, 1.0, d, n2, 1.0): hits += 1
    return hits
N = 6_000_000 if not QUICK else 1_500_000; box = 2.0; hits = mc_vex_disc(box, N, 11); p = hits / N
est = p * (2 * box)**3; se = (2 * box)**3 * math.sqrt(p * (1 - p) / N)
record("T4 thin-disc excluded volume = pi^2 R^3", abs(est - math.pi**2) < 4 * se, f"MC {est:.4f}+-{se:.4f} vs {math.pi**2:.4f}")

# ---- T5: cell-list completeness (the neighbour search in first_perc must cover every pair that can touch)
def cell_gaps(L, brad, pos):
    ncell = max(3, int(L // (2 * brad.max()))); cs = L / ncell
    cidx = (np.floor(pos / cs).astype(int)) % ncell
    n = len(pos); missed = 0; cand = 0
    for i in range(n):
        d = pos - pos[i]; d -= L * np.round(d / L); r2 = (d**2).sum(1); rr = brad + brad[i]
        c = np.where((r2 < rr**2) & (np.arange(n) > i))[0]; cand += len(c)
        dc = np.abs(cidx[c] - cidx[i]); dc = np.minimum(dc, ncell - dc)
        missed += int((dc.max(axis=1) > 1).sum()) if len(c) else 0
    return missed, cand, ncell, cs
def cell_gaps_big(L, R, hth, pos, ncell, cs):
    """Pairs that can touch (centre distance < sum of bounding radii) but lie in non-adjacent cells. Only pairs involving a
    large particle can be affected, so loop over the large ones."""
    brad = np.sqrt(R**2 + hth**2); cidx = (np.floor(pos / cs).astype(int)) % ncell; missed = 0; cand = 0
    big = np.where(brad > (cs - brad.max()) )[0]
    for i in big:
        d = pos - pos[i]; d -= L * np.round(d / L); r2 = (d**2).sum(1); rr = brad + brad[i]
        c = np.where(r2 < rr**2)[0]; c = c[c != i]; cand += len(c)
        dc = np.abs(cidx[c] - cidx[i]); dc = np.minimum(dc, ncell - dc); missed += int((dc.max(axis=1) > 1).sum())
    return missed // 1, cand
for lab, L, sig, nfac in (("monodisperse R=1 L=32", 32.0, 0.0, 2.5), ("lognormal sigma 0.3 (F12-like)", 32.0, 0.3, 2.5), ("lognormal sigma 0.4 (F4-like)", 32.0, 0.4, 2.5), ("lognormal sigma 0.5 (F5/F6-like)", 32.0, 0.5, 2.5)):
    n = int(nfac * L**3); R = np.ones(n) if sig == 0 else rng.lognormal(-0.5 * sig**2, sig, size=n)
    R = R[R <= 6.0]; hth = np.full(len(R), 0.02); pos = rng.random((len(R), 3)) * L
    brad = np.sqrt(R**2 + hth**2); ncell = max(3, int(L // (2 * brad.max()))); cs = L / ncell
    m, c = cell_gaps_big(L, R, hth, pos, ncell, cs)
    record("T5 cell list completeness: " + lab, m == 0, f"{m} candidate pairs missed of {c} examined for the large particles (ncell={ncell}, cell size {cs:.2f}, largest bounding radius {brad.max():.2f}, 2x = {2*brad.max():.2f})")

# ---- T6: percolation detector against an independent brute-force periodic-graph winding test (tolerance: identical Nc)
def contacts(mode, pos, rot, rad, hth, L):
    n = len(pos); E = []
    for i in range(n):
        d = pos - pos[i]; sh = -np.round(d / L); d = d + sh * L; r = np.sqrt((d**2).sum(1)); rr = np.sqrt(rad**2 + hth**2) + np.sqrt(rad[i]**2 + hth[i]**2)
        for j in np.where((r < rr) & (np.arange(n) > i))[0]:
            if mode == 0: hit = disc_contact(rot[i][:, 2], rad[i], d[j], rot[j][:, 2], rad[j])
            else: hit = gjk_cyl(np.zeros(3), rot[i][:, 2], rad[i], hth[i], d[j], rot[j][:, 2], rad[j], hth[j])
            if hit: E.append((i, j, tuple(sh[j].astype(int))))
    return E
def wraps(k, E, n):
    adj = [[] for _ in range(k)]
    for (i, j, s) in E:
        if i < k and j < k: adj[i].append((j, np.array(s))); adj[j].append((i, -np.array(s)))
    off = [None] * k
    for st in range(k):
        if off[st] is not None: continue
        off[st] = np.zeros(3, int); stack = [st]
        while stack:
            u = stack.pop()
            for (v, s) in adj[u]:
                w = off[u] + s
                if off[v] is None: off[v] = w; stack.append(v)
                elif (off[v] != w).any(): return True
    return False
def brute_Nc(E, n):
    lo, hi = 1, n
    if not wraps(n, E, n): return -1
    while lo < hi:
        mid = (lo + hi) // 2
        if wraps(mid, E, n): hi = mid
        else: lo = mid + 1
    return lo
mism = 0; tot = 0; det = []
for mode, A_, L, nreal in ((0, None, 6.0, 12 if not QUICK else 4), (1, 6.0, 6.0, 8 if not QUICK else 3)):
    for r in range(nreal):
        n = int(1.2 * L**3 * (1 if mode == 0 else 2)); pos = rng.random((n, 3)) * L; rot = random_rotations(n, rng)
        rad = np.ones(n); hth = np.zeros(n) if mode == 0 else np.full(n, 1.0 / A_)
        E = contacts(mode, pos, rot, rad, hth, L); a = run_once(mode, pos, rot, rad, hth, L, np.sqrt(rad**2 + hth**2)); b = brute_Nc(E, n)
        tot += 1; mism += (a != b)
        if a != b: det.append(f"mode{mode} real{r}: sim {a} vs brute {b}")
record("T6 first-wrapping-cluster detector vs brute force", mism == 0, f"{mism} mismatches in {tot} configurations " + "; ".join(det[:3]))

# ---- T7: orientation samplers and local-order definition
R3 = random_rotations(200000, rng)[:, :, 2]; c2 = (R3[:, 2]**2).mean(); P2 = (1.5 * R3[:, 2]**2 - 0.5).mean()
record("T7a isotropic orientations", abs(c2 - 1 / 3) < 0.003 and abs(P2) < 0.004, f"<cos^2>={c2:.4f} (1/3), <P2>={P2:+.4f} (0)")
nrm = runit(1000); Rm = rotation_with_normal(nrm, rng); orth = np.abs(np.einsum('nij,nik->njk', Rm, Rm) - np.eye(3)).max(); okn = np.abs(Rm[:, :, 2] - nrm).max()
record("T7b rotation_with_normal orthonormal, third column = normal", orth < 1e-12 and okn < 1e-12, f"max orthonormality error {orth:.1e}, normal error {okn:.1e}")
for kap, target in ((4.733, 0.5), (6.312, 0.6)):
    ct = fisher(kap, 400000, rng); S = (1.5 * ct**2 - 0.5).mean()
    record(f"T7c Fisher kappa={kap} gives S={target}", abs(S - target) < 0.004, f"S={S:.4f}")
# domain construction: scatter about each domain director must reproduce S when measured against the director
L = 36.0; ell = 6.0; npart = 60000; pos = rng.random((npart, 3)) * L; kap = 6.312
nseed = int(round(L**3 / ell**3)); rs = np.random.default_rng(5); seeds = rs.random((nseed, 3)) * L
nrm = aligned_normals(pos, L, ell, kap, np.random.default_rng(5))
# regenerate the identical director assignment to measure S against directors
rr = np.random.default_rng(5); seeds2 = rr.random((nseed, 3)) * L; dirs = rr.normal(size=(nseed, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
d = pos[:, None, :] - seeds2[None, :, :]; d -= L * np.round(d / L); idx = np.argmin((d**2).sum(2), axis=1)
Sdom = (1.5 * (nrm * dirs[idx]).sum(1)**2 - 0.5).mean(); Sbulk = np.linalg.eigvalsh(1.5 * (nrm[:, :, None] * nrm[:, None, :]).mean(0) - 0.5 * np.eye(3)).max()
record("T7d domain order parameter and isotropic bulk", abs(Sdom - 0.6) < 0.01 and Sbulk < 0.1, f"S about domain director {Sdom:.3f} (0.6), bulk order eigenvalue {Sbulk:.3f}")

# ---- T8: benchmark and finite-size scaling for thin discs (full result in fss_thin.json; quick check here)
if not QUICK:
    try:
        fs = json.load(open("fss_thin.json")); Ls = sorted(float(k) for k in fs)
        x = np.array([L**-1.141 for L in Ls if L >= 12]); y = np.array([np.mean(fs[str(L)]) for L in Ls if L >= 12]); w = np.array([1 / (np.std(fs[str(L)]) / math.sqrt(len(fs[str(L)])))**2 for L in Ls if L >= 12])
        Xm = np.vstack([np.ones_like(x), x]).T; W = np.diag(w); beta = np.linalg.solve(Xm.T @ W @ Xm, Xm.T @ W @ y); cov = np.linalg.inv(Xm.T @ W @ Xm)
        record("T8 thin-disc eta_c extrapolation (exponent 1/nu=1.141)", abs(beta[0] - 0.9614) < 3 * math.sqrt(cov[0, 0]) + 0.002, f"eta_inf = {beta[0]:.4f} +- {math.sqrt(cov[0,0]):.4f} vs 0.9614; L=32 ratio to benchmark {np.mean(fs['32.0'])/0.9614:.4f}")
    except Exception as e:
        record("T8 thin-disc extrapolation", False, "fss_thin.json missing: " + str(e))

json.dump(REP, open("test_report.json", "w"), indent=1)
print(f"\n{sum(r['passed'] for r in REP)} of {len(REP)} tests passed")
