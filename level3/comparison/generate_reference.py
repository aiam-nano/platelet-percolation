"""Generate the reference data a second implementation must reproduce. Run from the workdir (needs sim2.py, align.py).
Outputs (written next to this script): pairs.npz (deterministic contact vectors), configs.npz (deterministic first-percolation indices),
reference_stats.json (statistical reference values with standard errors)."""
import sys, os, math, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "code")))      # sim2.py, align.py
from sim2 import disc_contact, gjk_cyl, random_rotations, run_once, rotation_with_normal
from align import aligned_normals
# usage: python3 generate_reference.py [output_dir] [--no-stats]   (--no-stats skips the long statistical references)
args = [a for a in sys.argv[1:] if not a.startswith("--")]
out = os.path.abspath(args[0]) if args else HERE
os.makedirs(out, exist_ok=True)
# ---------------- deterministic pair vectors
rng = np.random.default_rng(20261024)
N = 4000
n1 = rng.normal(size=(N, 3)); n1 /= np.linalg.norm(n1, axis=1)[:, None]
n2 = rng.normal(size=(N, 3)); n2 /= np.linalg.norm(n2, axis=1)[:, None]
d = rng.uniform(-2.2, 2.2, size=(N, 3))
R1 = rng.uniform(0.5, 1.5, N); R2 = rng.uniform(0.5, 1.5, N)
h1 = rng.uniform(0.005, 0.6, N); h2 = rng.uniform(0.005, 0.6, N)      # half-thicknesses
z = np.zeros(3)
disc = np.array([disc_contact(n1[i], R1[i], d[i], n2[i], R2[i]) for i in range(N)])
cyl = np.array([gjk_cyl(z, n1[i], R1[i], h1[i], d[i], n2[i], R2[i], h2[i]) for i in range(N)])
np.savez(os.path.join(out, "pairs.npz"), n1=n1, n2=n2, d=d, R1=R1, R2=R2, h1=h1, h2=h2, disc=disc, cyl=cyl)
print("pairs: disc contacts", disc.sum(), "cylinder contacts", cyl.sum())
# ---------------- deterministic small configurations (box side 6, 120 particles, insertion order = array order)
cfg = {}
for k, (mode, t) in enumerate([(0, 0.0), (1, 0.25), (1, 0.05)]):
    idx = []
    for s in range(20):
        r = np.random.default_rng(1000*k+s); n = 120; L = 6.0
        pos = r.random((n, 3))*L; nrm = r.normal(size=(n, 3)); nrm /= np.linalg.norm(nrm, axis=1)[:, None]
        rot = rotation_with_normal(nrm, r); rad = np.ones(n); hth = np.full(n, t/2); brad = np.sqrt(rad**2+hth**2)
        Nc = run_once(mode, pos, rot, rad, hth, L, brad)
        cfg[f"m{mode}_t{t}_s{s}_pos"] = pos; cfg[f"m{mode}_t{t}_s{s}_nrm"] = nrm; idx.append(Nc)
    cfg[f"m{mode}_t{t}_Nc"] = np.array(idx); print("config set", mode, t, idx)
np.savez(os.path.join(out, "configs.npz"), **cfg)
# ---------------- statistical references (reduced quantities in the units of the spec)
def lognormal(n, s, rng, cap=6.0):
    R = rng.lognormal(-0.5*s*s, s, size=n); return R[R < cap]
def stat(case, L, nreal, seed):
    rng = np.random.default_rng(seed); vals = []
    for _ in range(nreal):
        n = int(2.5*L**3 if case != "disc" else 0.4*L**3)
        if case == "disc": R = np.ones(n); t = np.zeros(n); mode = 0
        elif case == "cyl_A8": R = np.ones(n); t = np.full(n, 2/8); mode = 1
        elif case == "cyl_A50": R = np.ones(n); t = np.full(n, 2/50); mode = 1
        elif case == "lognormal0.4_A80": R = lognormal(n, 0.4, rng); t = np.full(len(R), 2/80); mode = 1
        elif case == "stack3_A60": R = np.ones(n); t = np.full(n, 3*2/60); mode = 1
        elif case == "aligned_k6.312_ell6_A40": R = np.ones(n); t = np.full(n, 2/40); mode = 1
        n = len(R); pos = rng.random((n, 3))*L
        if case.startswith("aligned"):
            rot = rotation_with_normal(aligned_normals(pos, L, 6.0, 6.312, rng), rng)
        else: rot = random_rotations(n, rng)
        hth = t/2; brad = np.sqrt(R**2+hth**2)
        Nc = run_once(mode, pos, rot, R, hth, L, brad)
        if Nc > 0: vals.append((Nc/L**3*np.mean(math.pi*R[:Nc]**2*t[:Nc]) if mode == 1 else Nc/L**3*4*math.pi/3*np.mean(R[:Nc]**3)))
    v = np.array(vals); return dict(L=L, nreal=len(v), mean=float(v.mean()), se=float(v.std()/math.sqrt(len(v))), quantity=("n*V (nominal volume fraction)" if mode == 1 else "eta = (4pi/3) n R^3"))
ref = {}
if "--no-stats" in sys.argv: sys.exit(0)
for case in ["disc", "cyl_A8", "cyl_A50", "lognormal0.4_A80", "stack3_A60", "aligned_k6.312_ell6_A40"]:
    ref[case] = stat(case, 24.0, 100, 777+len(case)); print(case, ref[case], flush=True)
    json.dump(ref, open(os.path.join(out, "reference_stats.json"), "w"), indent=1)
