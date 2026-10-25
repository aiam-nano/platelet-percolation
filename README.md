# platelet-percolation

Monte Carlo simulator and simulation-calibrated law for the geometric percolation threshold of penetrable platelets of finite thickness.

**Status: private pre-release snapshot (v0.1.0, 25 October 2026, before the freeze).** Nothing here is licensed for reuse yet (see `LICENSE_PENDING.md`), the verification record is incomplete (see below), and the text in `docs/` consists of working documents carried over from the project that will be revised before v1.0.0. The planned route to a public release is described in `RELEASING.md`.

## What is in the repository

The simulator estimates the threshold of fully penetrable particles in a periodic cubic box. Particles are inserted one at a time at uniformly random positions and orientations, contacts are decided exactly, and clusters are tracked with a weighted union-find structure that carries periodic image offsets. The threshold of one realisation is the number of particles at which a cluster first closes on itself with a non-zero winding vector. Contact is exact for zero-thickness discs, finite-thickness cylinders (Gilbert–Johnson–Keerthi with analytic support functions), oblate and prolate spheroids (the same algorithm), and zero-thickness convex polygons (squares and equilateral triangles, by clipping each polygon to the plane of the other). Curved and rippled sheets use a mesh-edge-against-implicit-surface test, which is approximate: it misses 0.4 to 0.7 per cent of contacts at the production mesh spacing of one eighth of the radius (see `docs/15_Shape_Ripple_FSS_Checks_25Oct.md`).

The law, in `code/finallaw3.py` with its geometric constants in `code/mylaw.py`, gives the threshold volume fraction of finite-thickness platelets as a function of aspect ratio, restacking into tactoids, log-normal size dispersity, domain-correlated alignment and out-of-plane ripple. It reduces to the Yi–Tawerghi thin-disc value by construction. It has five fitted constants, which is more than the three the project set itself as a limit, and `code/calibrate_all.py` regenerates all of them from the stored simulation results in about one second.

```
code/                      simulator (sim2.py, sim3.py, sim4.py, sim2_axes.py), orientation generator (align.py),
                           the law (finallaw3.py, mylaw.py; finallaw2.py and finallaw.py kept as the earlier versions),
                           calibration (fitlaw.py, calibrate_all.py), exam generator (exam_sim.py)
code/verification/         regression suite (test_suite.py), finite-size-scaling runs, cross-validation, logs and results
code/literature_checks/    sim2e.py (adds the spheroid and polygon modes), polygon and spheroid contact checks,
                           thin-plate and open-face finite-size runs, ripple mesh test, experimental comparison
code/gap_tests/            the 24 pre-registered tests: predictions, drivers, scorer, raw results
code/harness/, code/exam/  scoring of the law on the sealed exam, and the exam truth with its SHA-256
docs/                      verification protocol, finite-size report, provenance ledger, exam scoring, shape and ripple checks
tests/check_reference.py   smoke test of the law and the calibration (also run by CI)
CHECKSUMS.sha256           SHA-256 of every file in this snapshot
```

## Quick start

Python 3.13 with the pinned packages in `requirements.txt` is the environment in which every result was produced (`ENVIRONMENT.txt`).

```
python3 -m pip install -r requirements.txt
python3 tests/check_reference.py                                   # law and calibration, a few seconds
cd code && PYTHONPATH=. python3 verification/test_suite.py quick  # 14 regression tests, about one minute
```

The first call of each compiled function takes 10 to 40 seconds because numba compiles it. Scripts import one another by module name, so run them from the folder that contains them with the parent `code/` folder on `PYTHONPATH`. `REPRODUCE.md` says how to regenerate each reported number.

## What has and has not been verified

The regression suite passes 14 of 15 tests; the one failure (T8) is the thin-disc threshold being 1.6 per cent above the published value, which is kept as a known discrepancy. The simulator's plateaus for zero-thickness discs, squares and triangles are 1.6, 1.7 and 2.9 per cent above the published thresholds of Yi and Tawerghi (2009), and its thin-limit spheroid threshold is about 30 per cent below a value that follows from Garboczi and co-workers' table; the causes are unknown. Every simulated threshold is divided by a plateau factor of 1.015 before it is used to fit or test the law, so the law matches the published thin-disc limit by construction. Of 24 pre-registered tests, 22 pass; the two failures (a bimodal radius mixture and a log-normal with width 0.8) are outside the law's stated scope. On a sealed exam the law's mean absolute error over the nine cases within its scope is 1.1 per cent, but the three aligned cases were mapped to the law's parameters after the truth was seen, so they are not blind.

The verification protocol has four levels. Levels 1 and 2 (the author's regression suite, and comparison with other groups' published numbers) have been run. Level 3 is an independent reimplementation by a separate AI system working only from a code-free specification, and Level 4 is an adversarial review by another AI system followed by human review. Level 3 has been run and its record will be added to the repository at the freeze; Level 4 and the human review are pending. A fresh sealed exam will be generated after the freeze.

## Use of AI

The code, the documents and most of the analysis were written by Claude (Anthropic) under the direction of the lead author, who defined the questions and the verification protocol. Because the person directing the work did not write the code line by line, the verification protocol includes independent reimplementation and adversarial review rather than relying on the author's own reading. The authorship list and the AI-use statement for the paper are to be settled by the authors.

## Citation and licence

`CITATION.cff` is a draft; the author list, version, date and DOI are to be decided. No licence has been chosen; see `LICENSE_PENDING.md`.
