"""Level 3 additions of 25 Oct 2026: deterministic inputs (written to ../reimplementation_cursor/inputs) and their reference answers
(written next to this script) for zero-thickness squares, zero-thickness equilateral triangles and oblate spheroids.
Needs sim2e.py and poly_test.py (literature_checks) on the path.   usage: python3 generate_reference_shapes.py <cursor_inputs_dir> [reference_output_dir]"""
import sys, os, json, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "code", "literature_checks")))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "code")))
from sim2e import poly_contact, gjk_ell, random_rotations, run_once, first_perc, EMPTY
from poly_test import polygon
inp = sys.argv[1]; out = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else HERE
N = 4000; z = np.zeros(3)
pairs_in, pairs_ref = {}, {}
HALF = {'square': 2.6, 'triangle': 3.2, 'spheroid': 2.2}
for k, shape in enumerate(('square', 'triangle', 'spheroid')):
    rng = np.random.default_rng(20261025 + k)
    F1 = random_rotations(N, rng); F2 = random_rotations(N, rng)              # columns: e1, e2, n
    d = rng.uniform(-HALF[shape], HALF[shape], size=(N, 3))
    rec = dict(e1_1=F1[:, :, 0], e2_1=F1[:, :, 1], n1=F1[:, :, 2], e1_2=F2[:, :, 0], e2_2=F2[:, :, 1], n2=F2[:, :, 2], d=d)
    if shape == 'spheroid':
        R1 = rng.uniform(0.5, 1.5, N); R2 = rng.uniform(0.5, 1.5, N)
        u1 = rng.choice([1.0, 3.0, 10.0, 100.0, 1000.0], N); u2 = rng.choice([1.0, 3.0, 10.0, 100.0, 1000.0], N)
        c1 = R1 / u1; c2 = R2 / u2; rec.update(R1=R1, R2=R2, c1=c1, c2=c2)
        hit = np.array([gjk_ell(z, F1[i, :, 2], R1[i], c1[i], d[i], F2[i, :, 2], R2[i], c2[i]) for i in range(N)])
    else:
        V = polygon(shape)
        hit = np.array([poly_contact(F1[i], d[i], F2[i], V) for i in range(N)])
    pairs_in[shape] = rec; pairs_ref[shape] = hit
    print(shape, 'contacts', int(hit.sum()), 'of', N)
np.savez(os.path.join(inp, 'shapes_pairs_inputs.npz'), **{f"{s}_{k}": v for s, r in pairs_in.items() for k, v in r.items()})
json.dump({s: {k: v.tolist() for k, v in r.items()} | {"n_pairs": N} for s, r in pairs_in.items()} | {"note": "for each shape, pair i: particle 1 at the origin with frame (e1_1, e2_1, n1) and particle 2 at d with frame (e1_2, e2_2, n2); the frames are orthonormal (columns of a rotation matrix). Polygons use e1, e2, n; the spheroid uses n only, with equatorial radii R1, R2 and polar semi-axes c1, c2."},
          open(os.path.join(inp, 'shapes_pairs_inputs.json'), 'w'))
np.savez(os.path.join(out, 'shapes_pairs_ref.npz'), **pairs_ref)
# configurations: box side 8, 250 particles, 20 per shape, insertion order = array order
cfg_in, cfg_ref = {}, {}
for k, (key, shape) in enumerate((('sq_L8', 'square'), ('tri_L8', 'triangle'), ('sph_L8', 'spheroid'))):
    lst, idx = [], []
    for s in range(20):
        r = np.random.default_rng(31000 + 1000 * k + s); n = 250; L = 8.0
        pos = r.random((n, 3)) * L; F = random_rotations(n, r)
        if shape == 'spheroid':
            rad = np.ones(n); hth = np.full(n, 1 / 1000.0); Nc = run_once(3, pos, F, rad, hth, L, np.ones(n) * 1.0000001)
        else:
            V = polygon(shape); br = float(np.linalg.norm(V, axis=1).max()); rad = np.ones(n); hth = np.zeros(n)
            Nc = first_perc(4, pos, F, rad, hth, L, max(3, int(L // (2 * br))), np.full(n, br), V, EMPTY['E'], 0, 0., 0., 0., 1., EMPTY['xtab'], EMPTY['utab'])
        lst.append(dict(pos=pos.tolist(), e1=F[:, :, 0].tolist(), e2=F[:, :, 1].tolist(), n=F[:, :, 2].tolist())); idx.append(int(Nc))
    cfg_in[key] = lst; cfg_ref[key] = np.array(idx); print(key, idx)
json.dump(cfg_in | {"note": "each configuration: 250 particles in a periodic box of side 8, in insertion order; pos, orthonormal frame (e1, e2, n); spheroid configurations use n only (R = 1, polar semi-axis 0.001); polygon configurations use the shape named by the key"},
          open(os.path.join(inp, 'shapes_configs_inputs.json'), 'w'))
np.savez(os.path.join(out, 'shapes_configs_ref.npz'), **cfg_ref)
np.savez(os.path.join(inp, 'shapes_configs_inputs.npz'), **{f"{k}_s{s}_{q}": np.array(c[q]) for k, lst in cfg_in.items() for s, c in enumerate(lst) for q in ('pos', 'e1', 'e2', 'n')})
