"""Smoke test of the law and its calibration. Run from the repository root:  python3 tests/check_reference.py

1. finallaw3.phi reproduces the reference values recorded at the 25 Oct 2026 snapshot (A = 10, 100, 1000 with every other factor off).
2. calibrate_all.py regenerates the five fitted constants from the stored simulation results.
Exit status 0 means both agree within the stated tolerances."""
import json, os, subprocess, sys
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
code = os.path.join(root, "code")
sys.path.insert(0, code)
os.chdir(code)
import finallaw3 as f

ok = True
for A, ref in ((10, 0.10865714), (100, 0.01396732), (1000, 0.00143743)):
    v = f.phi(A=A)
    good = abs(v / ref - 1) < 1e-5
    ok &= good
    print(f"{'PASS' if good else 'FAIL'} law A={A}: {v:.8f} (reference {ref})")

subprocess.run([sys.executable, "calibrate_all.py"], check=True, stdout=subprocess.DEVNULL)
res = json.load(open("calibration_result.json"))
expect = {"B_ (baseline)": [0.6628], "kappa (polydispersity)": [0.6248],
          "alignment (prefactor, gamma) free": [0.5878, 0.5763], "CRIP (ripple, thin-limit ratios)": [0.6474]}
for name, vals in expect.items():
    got = res[name]["value"]
    good = all(abs(g - e) < 5e-4 for g, e in zip(got, vals))
    ok &= good
    print(f"{'PASS' if good else 'FAIL'} {name}: {[round(g, 4) for g in got]} (reference {vals})")
sys.exit(0 if ok else 1)
