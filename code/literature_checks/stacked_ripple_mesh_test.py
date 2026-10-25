"""Independent test of the contact between two STACKED rippled sheets (thick platelets), action A4, 25 Oct 2026.

The production test (sim4.sheet_contact) treats a thick rippled platelet as M parametric sheets offset along the local z axis by
offs = linspace(-tv/2, tv/2, M) and reports contact when any layer of one crosses any layer of the other (a coarse edge mesh of one layer
against the implicit surface of the other).  This script repeats the test with nothing shared except the geometry definition: every layer
is a fine Delaunay triangle mesh (the construction of ripple_mesh_test.py), the layers are shifted along z, and two platelets touch when any
triangle of any layer of one meets any triangle of any layer of the other (exact segment-triangle tests).  It also runs a denser stack
(M_fine layers) to test whether M = 5 layers leave gaps through which a platelet can slip.

usage: python3 stacked_ripple_mesh_test.py npairs seed [case_index]"""
import sys, math, time, json
import numpy as np
import sim4
import ripple_mesh_test as rm

def stack(T, offs):
    return np.concatenate([T + np.array([0.0, 0.0, o])[None, None, :] for o in offs], axis=0)

def run(amp, lam, tv, M, M_fine, h_mesh, npairs, seed, sr=1.0, prod_h=(8, 16, 32)):
    rng = np.random.default_rng(seed)
    shs = {q: sim4.make_sheet('ripple', 1.0, h=1 / q, amp=amp, lam=lam) for q in prod_h}
    offs = np.linspace(-tv / 2, tv / 2, M); offf = np.linspace(-tv / 2, tv / 2, M_fine)
    T, X = rm.fine_mesh(amp, lam, h_mesh)
    brad = float(np.linalg.norm(X, axis=1).max()) + tv / 2
    rmax = brad * (1 + sr)
    S5 = stack(T, offs); Sf = stack(T, offf)
    cs = 4 * h_mesh + 0.05; lo = -(rmax + 0.2); nc = int(2 * (rmax + 0.2) / cs) + 1
    res = []
    for _ in range(npairs):
        while True:
            dv = (rng.random(3) - 0.5) * 2 * rmax
            if dv @ dv < rmax * rmax: break
        Ri, Rj = rm.sim_rand_rot(2, rng)
        hs = []
        for q in prod_h:
            sh = shs[q]; Y = np.empty((sh['V'].shape[0], 3)); g = np.empty(sh['V'].shape[0])
            hs.append(bool(sim4.sheet_contact(Ri, dv, Rj, sh['V'], sh['E'], sh['stype'], sh['sp0'], sh['sp1'], sh['sp2'], sh['R0'],
                                              sh['xtab'], sh['utab'], Y, g, offs, 1.0, sr)))
        A5 = np.empty_like(S5); B5 = np.empty_like(S5)
        rm.rotate_tris(S5, Ri, np.zeros(3), A5); rm.rotate_tris(S5 * sr if sr != 1.0 else S5, Rj, dv, B5)
        m5 = rm.mesh_contact(A5, B5, cs, lo, nc)
        Af = np.empty_like(Sf); Bf = np.empty_like(Sf)
        rm.rotate_tris(Sf, Ri, np.zeros(3), Af); rm.rotate_tris(Sf * sr if sr != 1.0 else Sf, Rj, dv, Bf)
        mf = rm.mesh_contact(Af, Bf, cs, lo, nc)
        res.append((hs, bool(m5), bool(mf)))
    return res, rmax

CASES = [  # name, amplitude, wavelength, total layer span tv, M, M_fine
    ("flat-like a0.05 thick tv0.2", 0.05 / (2 * math.pi) * 0 + 1e-9, 1.0, 0.2, 5, 11),
    ("slope2 lam1 tv0.2 (A=10)", 2.0 / (2 * math.pi), 1.0, 0.2, 5, 11),
    ("slope3 lam1 tv0.4 (A=5)", 3.0 / (2 * math.pi), 1.0, 0.4, 5, 11),
    ("slope3 lam0.7 tv0.4", 3.0 * 0.7 / (2 * math.pi), 0.7, 0.4, 5, 11),
    ("slope4 lam1.5 tv0.6 (A=3.3)", 4.0 * 1.5 / (2 * math.pi), 1.5, 0.6, 5, 11),
]

if __name__ == "__main__":
    npairs = int(sys.argv[1]); seed = int(sys.argv[2]); sel = [int(sys.argv[3])] if len(sys.argv) > 3 else range(len(CASES)); out = []
    for ci in sel:
        name, amp, lam, tv, M, Mf = CASES[ci]
        t0 = time.time(); res, rmax = run(amp, lam, tv, M, Mf, 1 / 20, npairs, seed)
        H = np.array([r[0] for r in res]); m5 = np.array([r[1] for r in res]); mf = np.array([r[2] for r in res]); vol = 4 / 3 * math.pi * rmax ** 3
        row = dict(case=name, tv=tv, M=M, M_fine=Mf, h_mesh=1 / 20, npairs=npairs, contacts_mesh_M=int(m5.sum()), contacts_mesh_Mfine=int(mf.sum()),
                   Vex_mesh_M=float(m5.mean() * vol), Vex_mesh_Mfine=float(mf.mean() * vol), Vex_err=float(vol * math.sqrt(m5.mean() * (1 - m5.mean()) / npairs)),
                   layers_gap_missed=int((mf & ~m5).sum()), secs=round(time.time() - t0))
        for c, q in enumerate((8, 16, 32)):
            row[f"prod_h{q}"] = dict(contacts=int(H[:, c].sum()), prod_only=int((H[:, c] & ~m5).sum()), mesh_only=int((~H[:, c] & m5).sum()),
                                     ratio_to_mesh=float(H[:, c].mean() / m5.mean()))
        print(json.dumps(row), flush=True); out.append(row)
    json.dump(out, open(f"stacked_ripple_mesh_test_{'_'.join(map(str, sel))}.json", "w"), indent=1)
