# Changelog

## v0.1.0, 25 October 2026 (private, pre-freeze snapshot)

First repository snapshot. Contents are the review package of 25 October 2026 (CHECKSUMS in the project folder record its 129 files) plus the repository files: `README.md`, `REPRODUCE.md`, `RELEASING.md`, `CHANGELOG.md`, `CITATION.cff`, `.zenodo.json`, `LICENSE_PENDING.md`, `ENVIRONMENT.txt`, `requirements.txt`, `.gitignore`, `.gitattributes`, `tests/check_reference.py`, a CI workflow and `CHECKSUMS.sha256`. `docs/01_Verification_Protocol.md` was added from the project folder.

Changes made to the files taken from the review package: session paths in the warning lines of `code/gap_tests/logB.txt` and `logC.txt` were removed; nothing else in the code, data or logs was altered. Generated outputs of the scripts (`calibration_result.json`, `fitpars.npy`, `code/test_report.json`) are excluded by `.gitignore`.

Checked on 25 October 2026 in the environment of `ENVIRONMENT.txt`: `tests/check_reference.py` passes, the quick regression suite passes 14 of 14 tests, and `regress_sim2e.py` finds no differences between `sim2.py` and `sim2e.py` in 36 configurations.

Not yet in the repository: the Level 3 record, the Level 4 findings, the human review, the fresh sealed exam, a licence, a DOI.
