# How to regenerate the numbers (v0.1.0, 25 Oct 2026, before the freeze)

Environment: `ENVIRONMENT.txt`. The first call of each numba function compiles it (10 to 40 seconds) and caches the result. The runs are single-threaded; two can run side by side on two cores. Every script takes its random seed from the command line or from a seed written in the script, so every result is deterministic given the seed, the code and the numba and numpy versions. Run scripts from inside the folder that contains them, with the `code/` folder on `PYTHONPATH` (for scripts in a sub-folder, `PYTHONPATH=..:.`). Check the files against `CHECKSUMS.sha256` with `sha256sum -c CHECKSUMS.sha256`.

## The law and its constants

`finallaw3.py` is the law (alignment prefactor 0.588; `finallaw2.py` with 0.625 is the record of the law frozen on 23 and 24 Oct, with which the sealed exam was first scored). `python3 calibrate_all.py` in `code/` regenerates every fitted constant from the stored simulation results in about one second and returns B 0.6628, kappa 0.6248, alignment prefactor 0.5878 with gamma 0.5763, and the ripple constant 0.6474 (chi-squared values in `code/verification/calibrate_all.log`). `python3 tests/check_reference.py` checks these and the law's reference values (A = 10, 100, 1000 with all other factors off: 0.10865714, 0.01396732, 0.00143743). `PYTHONPATH=. python3 verification/test_suite.py quick`, run in `code/`, runs the regression suite and takes about 40 seconds; it passes 14 of 14 tests. The full run adds the long finite-size test T8, which fails and is recorded as a known discrepancy (14 of 15 pass); it takes much longer and has not been repeated since 24 Oct (`code/verification/test_report.json`).

## Thin-disc offset and finite-size scaling

`verification/fss_thin.py` (nreal seed out.json) for boxes of 8 to 48 radii, `fss_thin_large.py` (L nreal seed out.json) for 64, 96 and 128, `axes_check.py` (L nreal seed out.json) for the single-axis criterion with `axes_fss_fit.py` for the fit, `fss_spheres.py` for the sphere control, `vex_check.py` for the excluded volume, `binary_check.py` for the binary mixture. The largest boxes take 10 to 40 minutes each. The stored results are the `.json` and `.log` files beside the scripts.

## Checks of 25 Oct (`code/literature_checks`)

`sim2e.py` is `sim2.py` with two added contact modes (3, oblate or prolate spheroid; 4, zero-thickness convex polygon); the other modes are unchanged. `python3 regress_sim2e.py` confirms the same first-wrap counts as `sim2.py` for the shared modes (36 configurations, no differences). `python3 poly_test.py` checks the polygon contact against the excluded volume A P / 2 (about 10 minutes). `./run_shapes.sh disc|square|triangle` runs `shape_fss.py` over boxes of 16 to 128 (about 30 minutes per shape on one core), `./run_span.sh disc` and `./run_span2.sh square|triangle` the open-boundary spanning runs, `./run_ell.sh` the spheroid finite-size runs, `python3 prolate_fss.py eps L nreal seed` the prolate runs, `python3 sphere_fss.py` the sphere control (about 4 minutes), and `python3 ripple_mesh_test.py 80000 20261025` the ripple mesh comparison (about 20 minutes). `python3 fss_analysis.py` summarises the jsonl files into `fss_analysis.txt`. `experimental_comparison.py` needs `finallaw3.py` and `mylaw.py` from `code/`.

## The pre-registered tests and the sealed exam

`code/gap_tests/gap_preds.py` wrote the predictions (`gap_preds.json`) before the simulations; the drivers `gap_cyl.py`, `gap_rip.py`, `gap_slab.py` and `gap_field.py` produce the `.pkl` results and `gap_score.py` scores them (22 of 24 inside their windows). `code/exam_sim.py` regenerates the exam truth from the seeds in the script; the truth file and its SHA-256 are in `code/exam/`, and `code/harness/score_finallaw_exam.py` scores the law against them. Neither needs simulation to re-score.

## Exam B, the offset diagnostics and the stacked-ripple test

`code/exam_B/predict_exam_B.py` writes the law's predictions for the ten cases of `spec.json` and the hashes of the law files (`predictions.json`, `predictions.sha256`); `exam_B_truth.py` simulates the truth (about 15 minutes on two cores) and writes a sealed file and its hash without printing values; `score_exam_B.py` checks every hash and scores. The truth file is not in the repository until the predictions have been scored; `criteria.md` says what counts as failure. `code/verification/offset_diagnostics.py` reruns the diagnostics of `docs/17_Offset_Diagnostics_25Oct.md` from the stored threshold files (`fss_*.json`, including the runs at 128 and 192 radii made with `fss_big.py`). `code/literature_checks/stacked_ripple_mesh_test.py 40000 7` repeats the independent mesh test for stacks of rippled sheets (about half an hour).

## What is not in this snapshot

The Co-Scientist scoring tools (see `RELEASING.md`). The Level 3 input arrays are now in `level3/`, and `tests/verify_level3_inputs.py` checks that the generator scripts reproduce them exactly. The simulations that produced the stored `.pkl`, `.json` and `.jsonl` files are long and were run once; the stored files are the record, and re-running should match them to within the statistical errors.
