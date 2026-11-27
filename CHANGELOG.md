# Changelog

## v1.0.0 (27 November 2026)

Experimental comparison. `code/literature_checks/experimental_comparison.py` and `.csv` now give the 16 cases from 11 papers with the transcribed values as used in Supplementary Table 4 of the paper (measured over model: median 4.4, range 0.14 to 177). No case is excluded. No code, model, calibration or simulation result changed, and the model files are byte-identical to v0.1.4; only the comparison script and its output differ.

Records. `code/exam_B/score_exam_B_output.txt` (the scoring of Exam B: ten of ten cases within 6 per cent, mean absolute error 1.6 per cent) is added to the repository.

Figures. `figures/` (new): the figures of the paper and the SI as PDF, the source data of every panel as CSV files, and the drawing scripts. Every figure regenerates from the repository alone, using the six small result files in `figures/data/` (`figures/README.md`).

Metadata. `CITATION.cff` and `.zenodo.json` name both authors and the version 1.0.0; `README.md` is updated. `CHECKSUMS.sha256` is regenerated. The repository is public at https://github.com/aiam-nano/platelet-percolation; the Zenodo archive has the concept DOI 10.5281/zenodo.24612913 (all versions) and, for this version, 10.5281/zenodo.24612914.

## v0.1.4 (private, pre-freeze)

Redaction. `tools/redact_transcripts.py --write` has been applied to the two session transcripts (`level3/cursor_instructions_execution_process.md`, `level4/cursor_instructions_execution_findings.md`): home-directory paths, the Cursor project-folder name and the operator's login name are replaced by placeholders (8 replacements). The unredacted originals are kept outside the repository. No code, law, calibration or result changed; `CHECKSUMS.sha256` is regenerated.

## v0.1.3 (private, pre-freeze)

Licence. `LICENSE` (MIT) and `NOTICE.md` (what the licence does not cover: outputs of Google's Hypothesis Generation tool and other AI services, extracted published values, third-party software) replace `LICENSE_PENDING.md`. `CITATION.cff` and `.zenodo.json` now name the licence. No code, law, calibration or result changed.

## v0.1.2, 25 October 2026, night (private, pre-freeze)

Diagnostics. `docs/17_Offset_Diagnostics_25Oct.md` and `code/verification/offset_diagnostics.py` (output `offset_diagnostics.out`): the estimator (Poisson ensemble, median), the crossing of wrapping curves, single-axis wrapping, tangent contacts and a free or fixed-limit finite-size fit all leave the thin-disc offset of 1.6 to 1.8 per cent in place. New data: 60 further realisations at 128 radii (0.9783 plus or minus 0.0010) and 40 at 192 radii (0.9778 plus or minus 0.0009), with `fss_big.py`; a slow approach to 0.9614 is excluded in practice.

Tests. `code/literature_checks/stacked_ripple_mesh_test.py` (new): independent triangle-mesh test of the contact between stacks of five rippled sheets, five cases of 40,000 pairs. At R/8 the production test finds 99.6 to 99.9 per cent of the mesh contacts, at R/32 99.7 to 100.1 per cent; eleven layers instead of five add 0.04 to 0.25 per cent. `tests/verify_level3_inputs.py` (new, in CI): regenerates the Level 3 input and reference arrays from the generator scripts and compares them bit for bit.

Records. `level3/package/inputs/` and `level3/comparison/*.npz` (the input and reference arrays, 5 MB in all) are now in the repository, and the generator scripts take an output directory. `tools/redact_transcripts.py` (new) removes local paths from the two session transcripts; it has not been applied.

Exam B. `code/exam_B/`: ten cases (four random, three aligned, three ripple ratios), predictions of `finallaw3.py` filed in `predictions.json` (SHA-256 `bc0b27f8c8268a6dd7f3e08f892073d515362cc01d80b9fe0f41f7e282047fd3`, 21:10 UTC), truth simulated afterwards and sealed (SHA-256 of the truth file `727b32b4ad4991aabd0520a0c210007d26c93a7f4898dda3872fc87e288cfea2` in `truth_sealed_B.sha256`; the file itself is not in the repository), scorer and criteria. The law, the calibration and the simulator are unchanged from v0.1.1: `tests/check_reference.py` passes with the same constants.

Still not in the repository: the human review, the Exam B truth file, a licence, a DOI.

## v0.1.1, 25 October 2026 (private, pre-freeze)

Changes after the Level 4 review (Run 2) and the Level 3 comparison (Run 1). Both records are now in the repository (`level3/`, `level4/`); the dispositions of the five Level 4 findings are in `docs/16_Level4_Dispositions_25Oct.md`.

Code. `code/literature_checks/sim2e.py`: (1) the sign test in `poly_interval` is now symmetric, so a polygon vertex or edge lying in the other polygon's plane gives the same answer from either side; (2) exactly parallel polygons are decided by a separating-axis test in the common plane (`poly_coplanar_overlap`) instead of always returning no contact; (3) `run_once` raises `ValueError` if the box side is less than four times the largest bounding radius, the condition under which the minimum-image rule is exact. `code/sim2.py`: the same guard in `run_once`. None of these changes any result for random orientations: 400,000 random square pairs gave identical decisions before and after, and threshold runs for squares, triangles and discs (six realisations each, box of 16 radii, seed 4) gave identical values.

Tests. `tests/test_polygon_contact.py` (new, in CI): 4,000 random coplanar pairs per shape against a polygon-clipping oracle, the vertex and edge cases from both sides and with roles swapped, and 3,000 random pairs per shape against an independent chord test. The regression suite (14 of 14), the larger-box test (T10), `regress_sim2e.py` (36 configurations, 0 differences) and `tests/check_reference.py` all pass with the changes.

Documents. Provenance ledger: the ripple factor Kell is described correctly (the inverse excluded-volume ratio of an elliptical footprint, not the arc-length ratio of a sinusoid). Document 02: the weighted disc mean covers the five boxes from L = 32 upward. New: `docs/16_Level4_Dispositions_25Oct.md`.

Records added. `level3/` (specification and operating instructions given to the second system, its code and report, the comparison script and reference generators, the comparison output) and `level4/` (findings, the reviewer's scripts, the reproduction outputs, the transcript).

Still not in the repository: the human review, the fresh sealed exam, a licence, a DOI, and the large Level 3 input arrays.

## v0.1.0, 25 October 2026 (private, pre-freeze snapshot)

First repository snapshot. Contents are the review package of 25 October 2026 (CHECKSUMS in the project folder record its 129 files) plus the repository files: `README.md`, `REPRODUCE.md`, `RELEASING.md`, `CHANGELOG.md`, `CITATION.cff`, `.zenodo.json`, `LICENSE_PENDING.md`, `ENVIRONMENT.txt`, `requirements.txt`, `.gitignore`, `.gitattributes`, `tests/check_reference.py`, a CI workflow and `CHECKSUMS.sha256`. `docs/01_Verification_Protocol.md` was added from the project folder.

Changes made to the files taken from the review package: session paths in the warning lines of `code/gap_tests/logB.txt` and `logC.txt` were removed; nothing else in the code, data or logs was altered. Generated outputs of the scripts (`calibration_result.json`, `fitpars.npy`, `code/test_report.json`) are excluded by `.gitignore`.

Checked on 25 October 2026 in the environment of `ENVIRONMENT.txt`: `tests/check_reference.py` passes, the quick regression suite passes 14 of 14 tests, and `regress_sim2e.py` finds no differences between `sim2.py` and `sim2e.py` in 36 configurations.

Not yet in the repository: the Level 3 record, the Level 4 findings, the human review, the fresh sealed exam, a licence, a DOI.
