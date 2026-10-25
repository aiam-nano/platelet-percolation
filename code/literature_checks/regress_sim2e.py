"""sim2e.py must reproduce sim2.py exactly for the modes they share (0 thin disc, 1 cylinder, 2 sheets): same configuration, same first-wrap count.
Run from the folder that holds both files (sim2.py from workdir)."""
import numpy as np, sim2, sim2e
bad = 0; n = 0
for seed in range(12):
    rng = np.random.default_rng(seed); L = 14.0; N = int(0.35 * L ** 3)
    pos = rng.random((N, 3)) * L; rot = sim2.random_rotations(N, rng); rad = np.ones(N)
    for mode, hth in ((0, np.zeros(N)), (1, np.full(N, 0.05))):
        br = np.sqrt(rad ** 2 + hth ** 2) + 1e-4
        a = sim2.run_once(mode, pos, rot, rad, hth, L, br); b = sim2e.run_once(mode, pos, rot, rad, hth, L, br)
        n += 1; bad += int(a != b)
    shS = sim2.make_sheet('ripple', 1.0, h=1 / 8, amp=0.2, lam=1.0); shE = sim2e.make_sheet('ripple', 1.0, h=1 / 8, amp=0.2, lam=1.0)
    br = np.full(N, shS['brad'])
    a = sim2.run_once(2, pos, rot, rad, np.zeros(N), L, br, sheet=shS); b = sim2e.run_once(2, pos, rot, rad, np.zeros(N), L, br, sheet=shE)
    n += 1; bad += int(a != b)
print(f"{n} configurations, {bad} differences between sim2.py and sim2e.py")
