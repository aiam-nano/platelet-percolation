"""Law, exam, polydispersity kernel, Kell, document arithmetic.
Run from review_work/: python test_law_docs.py
"""
import sys, os, json, math, pickle, hashlib
import numpy as np
from scipy.special import ellipe, ellipeinc

ROOT = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(ROOT, "code_copy")
sys.path.insert(0, CODE)
sys.path.insert(0, os.path.join(CODE, "literature_checks"))
sys.path.insert(0, os.path.join(CODE, "gap_tests"))
os.chdir(CODE)

import mylaw
import finallaw2 as F2
import finallaw3 as F3

def section(t):
    print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)

# ---------------------------------------------------------------- Kell vs arc-length ratio
section("Kell versus arc-length ratio")
def arc_ratio(a, n=200001):
    th = np.linspace(0, np.pi / 2, n)
    # uniform mean on [0, pi/2] equals (2/pi) * integral = period-average of sqrt(1+(a cos)^2)
    return float(np.mean(np.sqrt(1 + a * a * np.cos(th) ** 2)))

def kell_parts(a):
    th = np.linspace(0, np.pi / 2, 20001)
    k = 1 / np.mean(np.sqrt(1 + a * a * np.cos(th) ** 2))
    Kell = math.pi / (2 * k * ellipe(1 - k * k))
    # closed form of the arc ratio
    # (2/pi) * sqrt(1+a^2) * E( a^2/(1+a^2) )
    closed = (2 / math.pi) * math.sqrt(1 + a * a) * ellipe(a * a / (1 + a * a))
    return k, Kell, closed

for a in (0.0, 0.5, 1.0, 1.2, 1.5, 2.0, 3.0):
    if a == 0:
        print(f"a={a:.1f}  arc={arc_ratio(0.0):.6f}  Kell(limit) n/a")
        continue
    k, Kell, closed = kell_parts(a)
    ar = arc_ratio(a)
    print(f"a={a:.1f}  arc_quad={ar:.6f}  arc_ellip={closed:.6f}  Kell={Kell:.6f}  Kell/arc={Kell/ar:.6f}  k={k:.6f}")

# ---------------------------------------------------------------- polydispersity kernel
section("polydispersity matrix versus Isihara kernel")
pi = math.pi

def Vex(R1, T1, R2, T2):
    """Orientation-averaged excluded volume of two cylinders. T is full thickness."""
    V1 = pi * R1 ** 2 * T1; V2 = pi * R2 ** 2 * T2
    S1 = 2 * pi * R1 ** 2 + 2 * pi * R1 * T1; S2 = 2 * pi * R2 ** 2 + 2 * pi * R2 * T2
    M1 = pi * T1 + pi ** 2 * R1; M2 = pi * T2 + pi ** 2 * R2
    return V1 + V2 + (S1 * M2 + S2 * M1) / (4 * pi)

def moments_ln(sig, k):
    return math.exp(k * (k - 1) * sig * sig / 2.0)

def code_lambda(sig, A):
    return mylaw.lam_max(sig, A)

def correct_lambda(sig, T, ngrid=96):
    """Largest eigenvalue of the excluded-volume operator, number-mean radius 1, fixed full thickness T.
    Returns (lambda, mean particle volume)."""
    if sig < 1e-8:
        R = np.array([1.0]); w = np.array([1.0])
    else:
        # Gauss-Hermite: int f(x) e^{-x^2} dx.  Standard normal u = x*sqrt(2).
        x, wh = np.polynomial.hermite.hermgauss(ngrid)
        u = x * math.sqrt(2.0)
        w = wh / math.sqrt(math.pi)   # integrates to 1 against the standard normal
        mu = -0.5 * sig * sig
        R = np.exp(mu + sig * u)
    # matrix K_ij = Vex(Ri,T,Rj,T) * w_j   eigenvalue of K
    RR1 = R[:, None]; RR2 = R[None, :]
    K = Vex(RR1, T, RR2, T) * w[None, :]
    lam = max(np.linalg.eigvals(K).real)
    Vp = np.sum(w * pi * R ** 2 * T)
    return lam, Vp

def report_poly(sig, A):
    T_full = 2.0 / A          # A = diameter/thickness, mean radius 1
    T_code = 1.0 / A          # what mylaw.lam_max uses internally
    lm = code_lambda(0.0, A); ls = code_lambda(sig, A)
    # code matrix factor, without the empirical (1+kappa sig^2)
    # P = (lm * exp(sig^2) / ls) * (1+kappa sig^2)
    mfac = (lm * math.exp(sig * sig) / ls) if sig > 0 else 1.0
    lam0, Vp0 = correct_lambda(0.0, T_full)
    lam, Vp = correct_lambda(sig, T_full)
    # monodisperse reference at R=1
    correct_fac = (Vp / lam) / (Vp0 / lam0)
    lam0b, Vp0b = correct_lambda(0.0, T_code)
    lamb, Vpb = correct_lambda(sig, T_code)
    fac_halfT = (Vpb / lamb) / (Vp0b / lam0b)
    Pcode = mylaw.P(sig, A, 0.625)
    print(f"sig={sig:.2f} A={A:g}: code_matrix_factor={mfac:.4f}  P(kappa=0.625)={Pcode:.4f}  "
          f"Isihara T=2/A factor={correct_fac:.4f}  Isihara T=1/A factor={fac_halfT:.4f}  "
          f"code/Isihara(T=2/A)={mfac/correct_fac:.4f}")

print("monodisperse code lambda vs Vex scale check at R=1")
for A, Tname in ((1e6, "thin"), (100, "A100"), (30, "A30"), (10, "A10"), (8, "A8")):
    T = 2 / A
    vex = Vex(1, T, 1, T)
    lc = code_lambda(0.0, A)
    # sum of code coefficients
    print(f"  {Tname}: Vex={vex:.6f}  code_lam={lc:.6f}  Vex/code_lam={vex/lc:.3f}")

print("\nfactor comparison (1 = code matrix matches Isihara ratio)")
for sig in (0.2, 0.3, 0.4, 0.5, 0.8):
    for A in (1e6, 80, 60, 30, 20, 10, 8):
        report_poly(sig, A)

# calibration thin ratios versus full P
section("calibration targets versus P()")
# from fitlaw.pol and fitlaw.fin
pol = [(1e6, .1, 1.007), (1e6, .2, 1.021), (1e6, .3, 1.053), (1e6, .4, 1.099), (1e6, .45, 1.136)]
fin = [(33.333, .3, 0.843), (30, .4, 0.797), (10, .3, 0.925)]
print("thin Br vs (1+kappa s^2) vs full P")
for A, s, Br in pol:
    drift = 1 + 0.625 * s * s
    print(f"  s={s:.2f} Br_sim={Br:.3f} drift={drift:.3f} ({100*(drift/Br-1):+.1f}%)  fullP={mylaw.P(s,A,0.625):.3f} ({100*(mylaw.P(s,A,0.625)/Br-1):+.1f}%)")
print("finite-aspect sim ratio vs full P")
for A, s, y in fin:
    print(f"  A={A:g} s={s:.2f} sim={y:.3f} P={mylaw.P(s,A,0.625):.3f} ({100*(mylaw.P(s,A,0.625)/y-1):+.1f}%)")

# ---------------------------------------------------------------- lognormal cap
section("lognormal truncation (mean radius, cap, no redraw)")
def truncated_stats(sig, rcap, n=2_000_000, seed=1):
    rng = np.random.default_rng(seed)
    R = rng.lognormal(-0.5 * sig * sig, sig, size=n)
    frac = np.mean(R > rcap)
    Rc = R[R <= rcap]
    return frac, Rc.mean(), np.mean(Rc ** 2), np.mean(Rc ** 3)

for sig, cap in ((0.3, 6), (0.4, 6), (0.5, 6), (0.8, 5), (0.8, 6)):
    frac, m, m2, m3 = truncated_stats(sig, cap)
    # ideal moments
    print(f"sig={sig} cap={cap}: P(R>cap)={frac:.6f}  E[R|cap]={m:.6f} (ideal 1)  "
          f"E[R^2|cap]={m2:.4f} (ideal {math.exp(sig*sig):.4f})  "
          f"volume bias E[R^2|cap]/E[R^2]={m2/math.exp(sig*sig):.4f}")

# ---------------------------------------------------------------- exam scores
section("exam scores finallaw2 and finallaw3")
raw = open(os.path.join(CODE, "exam", "truth_sealed.json"), "rb").read()
print("sha256", hashlib.sha256(raw).hexdigest())
t = json.loads(raw)
# Use the generator measurement if align is importable; also score with the documented S values
S_doc = {0.5: 0.1662, 0.6: 0.2394}
def score(F, Sloc, label):
    cases = {
        "F1": dict(A=50), "F2": dict(A=8), "F3": dict(A=60, N=4), "F4": dict(A=80, sig=.4),
        "F8": dict(A=40, S=Sloc[0.6], d=6), "F9": dict(A=60, N=3, S=Sloc[0.6], d=6),
        "F12": dict(A=60, N=3, sig=.3, S=Sloc[0.5], d=6),
    }
    errs = []
    for k, a in cases.items():
        occ = 1 - math.exp(-F.phi(**a)); e = 100 * (occ / t[k] - 1); errs.append(abs(e))
        print(f"  {label} {k:4s} law {occ:.5f} truth {t[k]:.5f} err {e:+.2f}%")
    for k in ("F10", "F11"):
        r = F.Krip(1.2, 1.0); e = 100 * (r / t[k] - 1); errs.append(abs(e))
        print(f"  {label} {k:4s} ratio {r:.4f} truth {t[k]:.4f} err {e:+.2f}%")
    print(f"  {label} MAE {np.mean(errs):.3f}% over {len(errs)}")
    return errs

print("with documented S_loc")
score(F2, S_doc, "law2")
score(F3, S_doc, "law3")

# limits
section("law limits and pole")
for name, F in (("law2", F2), ("law3", F3)):
    # thin: eta = phi * 2A/3
    for A in (1e3, 1e6):
        ph = F.phi(A=A); eta = ph * 2 * A / 3
        print(f"  {name} A={A:g} phi={ph:.6e} eta={eta:.6f}")
    ph = F.phi(A=1e6, sig=0, S=0)
    print(f"  {name} monodisperse zero-order thin eta {ph*2e6/3:.6f}")
# pole of Kal
print("finallaw3 Kal pole scan")
for S, d in ((0.5, 2), (1.0, 2), (0.5, 1), (1.0, 1), (0.8, 1.5)):
    try:
        k = F3.Kal(S, d)
        print(f"  S={S} d={d} Kal={k:.4f}")
    except Exception as e:
        print(f"  S={S} d={d} EXC {e}")

# Ae stack check
section("stack Ae versus single cylinder")
# law phi(A=60,N=4) should equal phi(A=15) for the base (no poly, no align)
print("phi(A=60,N=4)", F3.phi(A=60, N=4), "phi(A=15)", F3.phi(A=15))
print("phi(A=60,N=3)", F3.phi(A=60, N=3), "phi(A=20)", F3.phi(A=20))

# ---------------------------------------------------------------- experimental ratios
section("experimental comparison recompute")
csv_path = os.path.join(CODE, "literature_checks", "experimental_comparison.csv")
import csv
rows = list(csv.DictReader(open(csv_path)))
n_ok = 0
mids = []
for r in rows:
    Alo, Ahi = float(r["A_low"]), float(r["A_high"]); sig = float(r["sigma_lnR"])
    tlo, thi = float(r["measured_low_vol%"]), float(r["measured_high_vol%"])
    Am = math.sqrt(Alo * Ahi)
    pm = 100 * F3.phi(Am, sig=sig)
    plo = 100 * F3.phi(Ahi, sig=sig); phi_loA = 100 * F3.phi(Alo, sig=sig)
    tm = math.sqrt(tlo * thi)
    ratio = tm / pm
    rlo = tlo / phi_loA; rhi = thi / plo
    ok = (abs(ratio - float(r["ratio_mid"])) < 0.02 and abs(float(r["law_vol%_at_A_high"]) - plo) < 0.002)
    n_ok += ok
    mids.append((r["system"][:40], ratio, float(r["ratio_mid"]), ok))
print(f"rows matching csv mid-ratio and law_at_A_high: {n_ok}/{len(rows)}")
for s, ratio, csv_r, ok in mids:
    flag = "" if ok else "  MISMATCH"
    print(f"  {s:40s} recomputed {ratio:.2f} csv {csv_r:.2f}{flag}")
between = [m for m in mids if 0.7 <= m[1] <= 1.3]
print(f"mid-ratios in [0.7,1.3]: {len(between)} -> {[b[0] for b in between]}")
print(f"min mid {min(m[1] for m in mids):.2f}  max mid {max(m[1] for m in mids):.2f}")

# ---------------------------------------------------------------- gap score (imports from gap_tests cwd)
section("gap_score")
os.chdir(os.path.join(CODE, "gap_tests"))
import gap_score  # prints

# ---------------------------------------------------------------- document plateaus and median extrapolation
section("document numbers from raw jsonl / json")
os.chdir(os.path.join(CODE, "literature_checks"))

def wfit(x, y, e):
    A = np.vander(x, 2, increasing=True); W = np.diag(1 / np.array(e) ** 2)
    cov = np.linalg.inv(A.T @ W @ A); b = cov @ A.T @ W @ y
    chi2 = float(((A @ b - y) / e) @ ((A @ b - y) / e))
    return b, np.sqrt(np.diag(cov)), chi2

span = [json.loads(l) for l in open("span_open.jsonl")]
for sh in ("disc", "square", "triangle"):
    r = sorted([x for x in span if x["shape"] == sh], key=lambda x: x["L"])
    L = np.array([x["L"] for x in r], float)
    for key in ("eta_axis_mean", "eta_median"):
        y = np.array([x[key] for x in r]); e = np.array([x["eta_axis_err"] for x in r])
        # median has no separate stderr in the file; reuse axis err as the author's likely choice, and also an unweighted fit
        b, s, chi = wfit(1 / L, y, e)
        print(f"  {sh} linear 1/L of {key}: intercept {b[0]:.4f} +- {s[0]:.4f} chi2 {chi:.2f}/{len(r)-2}")
    print("   L, mean, median:", [(int(x["L"]), round(x["eta_axis_mean"], 4), round(x["eta_median"], 4), x["n"]) for x in r])

# sphere log recompute
section("sphere fss weighted means")
sph = [(10, 0.3553, 0.0037), (14, 0.3412, 0.0042), (20, 0.3442, 0.0031), (28, 0.3450, 0.0027), (40, 0.3397, 0.0023), (56, 0.3417, 0.0022)]
for Lmin in (10, 14, 20):
    sub = [t for t in sph if t[0] >= Lmin]
    w = np.array([1 / e ** 2 for _, _, e in sub]); y = np.array([m for _, m, _ in sub])
    print(f"  weighted mean L>={Lmin}: {np.sum(w*y)/np.sum(w):.4f} +- {1/math.sqrt(np.sum(w)):.4f}")
L = np.array([t[0] for t in sph], float); y = np.array([t[1] for t in sph]); e = np.array([t[2] for t in sph])
b, s, chi = wfit(L ** (-1 / 0.8765), y, e)
print(f"  fit all L: {b[0]:.4f} +- {s[0]:.4f} chi2 {chi:.2f}/{len(sph)-2}")

# old fss thin from verification json
section("doc 02 fss_thin.json check")
os.chdir(os.path.join(CODE, "verification"))
fs = json.load(open("fss_thin.json"))
def mean_se(key):
    a = np.array(fs[key], float)
    return a.mean(), a.std(ddof=1) / math.sqrt(len(a)), len(a)
for k in sorted(fs, key=lambda s: float(s)):
    m, se, n = mean_se(k)
    print(f"  L={k}: mean {m:.4f} +- {se:.4f} n={n}")
keys = [k for k in fs if float(k) >= 32]
# four largest and five largest
for label, ks in (("L>=32", keys), ("four largest", sorted(keys, key=float)[-4:]), ("five largest", sorted(keys, key=float)[-5:])):
    w = []; y = []
    for k in ks:
        m, se, n = mean_se(k); w.append(1 / se ** 2); y.append(m)
    w = np.array(w); y = np.array(y)
    print(f"  weighted {label} {ks}: {np.sum(w*y)/np.sum(w):.4f} +- {1/math.sqrt(w.sum()):.4f}")

# ripple ratios at production h/8 versus the fine mesh (h=0.025)
section("ripple mesh ratios")
os.chdir(os.path.join(CODE, "literature_checks"))
data = json.load(open("ripple_mesh_test.json"))
for d in data:
    if abs(d["h_mesh"] - 0.025) > 1e-9:
        continue
    h8 = d["prod_h8"]["ratio_to_mesh"]; h32 = d["prod_h32"]["ratio_to_mesh"]
    print(f"  {d['case']:22s} fine-mesh contacts {d['contacts_mesh']:6d}  prod/8 {h8:.4f} ({100*(h8-1):+.2f}%)  prod/32 {h32:.4f} ({100*(h32-1):+.2f}%)  "
          f"mesh_only h8 {d['prod_h8']['mesh_only']} prod_only {d['prod_h8']['prod_only']}")

print("\nDONE")
