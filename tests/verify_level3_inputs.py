"""Check that the Level 3 inputs and reference answers stored in level3/ are regenerated exactly by the generator scripts.
Runs level3/comparison/generate_reference.py (deterministic part only) and generate_reference_shapes.py into a temporary folder and compares
every array with the stored files, bit for bit.  The statistical references (reference_stats.json) are long simulations and are not rerun.
Run from the repository root:  python tests/verify_level3_inputs.py      (a few seconds; exit status 1 on any difference)"""
import os, subprocess, sys, tempfile, numpy as np
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cmp_dir = os.path.join(root, "level3", "comparison"); inp_dir = os.path.join(root, "level3", "package", "inputs")
bad = 0
def same(a, b, label):
    global bad
    ok = a.shape == b.shape and np.array_equal(a, b); bad += (not ok)
    return ok
with tempfile.TemporaryDirectory() as tmp:
    ref_out = os.path.join(tmp, "ref"); shp_in = os.path.join(tmp, "inputs"); os.makedirs(shp_in)
    subprocess.run([sys.executable, os.path.join(cmp_dir, "generate_reference.py"), ref_out, "--no-stats"], check=True, stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable, os.path.join(cmp_dir, "generate_reference_shapes.py"), shp_in, ref_out], check=True, stdout=subprocess.DEVNULL)
    def cmp_files(new_path, old_path, label, only=None, skip=()):
        n, o = np.load(new_path), np.load(old_path); keys = [k for k in o.keys() if k not in skip]
        miss = [k for k in keys if k not in n.keys()]
        diff = [k for k in keys if k in n.keys() and not same(n[k], o[k], k)]
        extra = [k for k in n.keys() if k not in o.keys() and k not in skip]
        status = "IDENTICAL" if not (miss or diff) else f"DIFFERENT (missing {len(miss)}, differing {len(diff)})"
        print(f"  {label}: {len(keys)} arrays compared, {status}" + (f"; generated file also holds {len(extra)} arrays not in the stored file (answers withheld from the package)" if extra else ""))
        return not (miss or diff)
    ok = True
    ok &= cmp_files(os.path.join(ref_out, "pairs.npz"), os.path.join(cmp_dir, "pairs.npz"), "comparison/pairs.npz (inputs and reference answers)")
    ok &= cmp_files(os.path.join(ref_out, "pairs.npz"), os.path.join(inp_dir, "pairs_inputs.npz"), "package/inputs/pairs_inputs.npz (inputs only)")
    ok &= cmp_files(os.path.join(ref_out, "configs.npz"), os.path.join(cmp_dir, "configs.npz"), "comparison/configs.npz")
    ok &= cmp_files(os.path.join(ref_out, "configs.npz"), os.path.join(inp_dir, "configs_inputs.npz"), "package/inputs/configs_inputs.npz (inputs only)")
    ok &= cmp_files(os.path.join(shp_in, "shapes_pairs_inputs.npz"), os.path.join(inp_dir, "shapes_pairs_inputs.npz"), "package/inputs/shapes_pairs_inputs.npz")
    ok &= cmp_files(os.path.join(shp_in, "shapes_configs_inputs.npz"), os.path.join(inp_dir, "shapes_configs_inputs.npz"), "package/inputs/shapes_configs_inputs.npz")
    ok &= cmp_files(os.path.join(ref_out, "shapes_pairs_ref.npz"), os.path.join(cmp_dir, "shapes_pairs_ref.npz"), "comparison/shapes_pairs_ref.npz")
    ok &= cmp_files(os.path.join(ref_out, "shapes_configs_ref.npz"), os.path.join(cmp_dir, "shapes_configs_ref.npz"), "comparison/shapes_configs_ref.npz")
print("ALL LEVEL 3 INPUTS AND REFERENCES REGENERATE EXACTLY" if ok else "MISMATCH")
sys.exit(0 if ok else 1)
