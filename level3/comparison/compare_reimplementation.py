"""Compare the independent reimplementation (results.json from the Cursor folder) with the reference simulator.
Usage: python3 compare_reimplementation.py path/to/results.json
Run it from anywhere; it reads the reference files in this folder (pairs.npz, configs.npz, reference_stats.json, reference_axes/).
Sections: A pairs; B configs; C statistics (pass/fail, 3 combined standard errors + systematic allowance); D internal consistency of the raw
arrays; E diagnostics (per-axis winding and alternative estimators; NOT pass/fail, for interpreting the 1.6 per cent thin-disc offset)."""
import sys, json, os, math, numpy as np
here = os.path.dirname(os.path.abspath(__file__)); r = json.load(open(sys.argv[1]))
P = np.load(os.path.join(here, "pairs.npz")); C = np.load(os.path.join(here, "configs.npz")); ref = json.load(open(os.path.join(here, "reference_stats.json")))
ok = True; K = 4*math.pi/3
print("== A. Deterministic pair tests =="); 
for k in ("disc", "cyl"):
    mism = int((np.array(r["pairs"][k], bool) != P[k]).sum()); lim = int(0.002*len(P[k])); good = mism <= lim; ok &= good
    print(f"pairs {k:4s}: {mism} mismatches of {len(P[k])} (allowed {lim}) {'PASS' if good else 'FAIL'}")
print("== B. Deterministic configuration tests ==")
for key, v in ((k, r["configs"][k]) for k in ("m0_t0.0", "m1_t0.25", "m1_t0.05") if k in r["configs"]):
    mism = int((np.array(v) != C[key+"_Nc"]).sum()); good = mism <= 1; ok &= good
    print(f"config {key:10s}: {mism} of 20 differ (allowed 1) {'PASS' if good else 'FAIL'}")
print("== C. Statistics ==")
def cmp(label, m, e, s):
    global ok
    tol = 3*math.hypot(e, s["se"]) + s.get("sys", 0.005)*s["mean"]; dif = m - s["mean"]; good = abs(dif) <= tol; ok &= good
    print(f"{label:30s}: ref {s['mean']:.5f}+-{s['se']:.5f}  new {m:.5f}+-{e:.5f}  diff {100*dif/s['mean']:+.2f}%  tol {100*tol/s['mean']:.2f}%  {'PASS' if good else 'FAIL'}")
St = r["stats"]
mapping = {"disc": "disc", "disc_L32": "disc_L32_highprecision", "cyl_A8": "cyl_A8", "cyl_A50": "cyl_A50", "stack3_A60": "stack3_A60", "lognormal0.4_A80": "lognormal0.4_A80", "aligned_k6.312_ell6_A40": "aligned_k6.312_ell6_A40"}
for case, refkey in mapping.items():
    if case not in St: print(f"{case}: MISSING"); ok = False; continue
    cmp(case, St[case]["mean"], St[case]["se"], ref[refkey])
if "disc_L48" in St and "disc_L64" in St:
    w = [1/St[c]["se"]**2 for c in ("disc_L48", "disc_L64")]; m = sum(wi*St[c]["mean"] for wi, c in zip(w, ("disc_L48", "disc_L64")))/sum(w); e = 1/math.sqrt(sum(w))
    cmp("disc_L48+L64 pooled", m, e, ref["disc_plateau_L48_L64_pooled"])
else: print("disc_L48 / disc_L64: MISSING"); ok = False
print("== A2/B2/C2. Shape additions of 25 Oct (square, triangle, spheroid) ==")
SP = np.load(os.path.join(here, "shapes_pairs_ref.npz")); SC = np.load(os.path.join(here, "shapes_configs_ref.npz"))
for k in ("square", "triangle", "spheroid"):
    if k not in r["pairs"]: print(f"pairs {k}: MISSING"); ok = False; continue
    mism = int((np.array(r["pairs"][k], bool) != SP[k]).sum()); lim = int(0.002*len(SP[k])); good = mism <= lim; ok &= good
    print(f"pairs {k:8s}: {mism} mismatches of {len(SP[k])} (allowed {lim}) {'PASS' if good else 'FAIL'}")
for key in ("sq_L8", "tri_L8", "sph_L8"):
    if key not in r["configs"]: print(f"config {key}: MISSING"); ok = False; continue
    mism = int((np.array(r["configs"][key]) != SC[key]).sum()); good = mism <= 1; ok &= good
    print(f"config {key:8s}: {mism} of 20 differ (allowed 1) {'PASS' if good else 'FAIL'}")
for case in ("square_L32", "square_L48", "triangle_L32", "triangle_L48", "spheroid1000_L32", "spheroid1000_L48"):
    if case not in St: print(f"{case}: MISSING"); ok = False; continue
    cmp(case, St[case]["mean"], St[case]["se"], ref[case])
print("== D. Internal consistency of raw arrays ==")
Ls = {"disc": 24, "disc_L32": 32, "disc_L48": 48, "disc_L64": 64, "square_L32": 32, "square_L48": 48, "triangle_L32": 32, "triangle_L48": 48}
for c, L in Ls.items():
    raw = r.get("raw", {}).get(c)
    if not raw: print(f"{c}: raw arrays MISSING"); ok = False; continue
    a = {k: np.array(raw[k], float) for k in ("N_any", "N_x", "N_y", "N_z")}
    viol = int((a["N_any"] != np.minimum(np.minimum(a["N_x"], a["N_y"]), a["N_z"])).sum())
    eta = K*a["N_any"].mean()/L**3; se_raw = K*a["N_any"].std(ddof=1)/math.sqrt(len(a["N_any"]))/L**3
    stat_ok = abs(eta - St[c]["mean"]) < 1e-6*max(1, abs(eta)) + 1e-9; ok &= (viol == 0 and stat_ok)
    print(f"{c:9s}: N_any != min(N_x,N_y,N_z) in {viol} realisations; mean eta from raw {eta:.5f} (+-{se_raw:.5f}) vs stats {St[c]['mean']:.5f} -> {'consistent' if stat_ok and viol == 0 else 'INCONSISTENT'}")
print("== E. Diagnostics (no pass/fail): per-axis winding, reduced density = (4pi/3) N / L^3 ==")
for c, L in Ls.items():
    raw = r.get("raw", {}).get(c)
    if not raw: continue
    line = f"{c:9s} (L={L}) new :"
    for k in ("N_any", "N_x", "N_y", "N_z"):
        x = K*np.array(raw[k], float)/L**3; line += f"  {k[2:]} {x.mean():.4f}+-{x.std(ddof=1)/math.sqrt(len(x)):.4f}"
    print(line)
    f = os.path.join(here, "reference_axes", f"axes_L{L}.json")
    if c.startswith("disc") and os.path.exists(f):
        a = np.array(json.load(open(f)), float)*K/L**3; line = f"{'':9s}  ref :"
        for j, k in enumerate(("any", "x", "y", "z")): line += f"  {k} {a[:,j].mean():.4f}+-{a[:,j].std(ddof=1)/math.sqrt(len(a)):.4f}"
        print(line)
print("OVERALL", "PASS" if ok else "FAIL (read README_COMPARISON.md for the decision logic)")
