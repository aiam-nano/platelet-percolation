# Route to v1.0.0

Working note for the authors. Nothing here has been done unless it is marked done.

## Where things stand

v0.1.0 (25 October 2026) was a private, pre-freeze snapshot; v0.1.1 (the same day) adds the Level 3 and Level 4 records and the polygon-contact correction; v0.1.2 (the same night) adds the offset diagnostics, the stacked-ripple mesh test, the Level 3 input arrays and Exam B. It is a private, pre-freeze snapshot of the simulator, the law, the calibration, the regression suite, the verification scripts and the reproduction note. It has no licence and no DOI.

## Before the freeze (v0.x)

Level 4 is settled (see `docs/16_Level4_Dispositions_25Oct.md`); the independent human review is still to be done. Decide how to state the 0.4 to 0.7 per cent ripple-contact bias and whether to refit the ripple constant at a finer mesh. The layered ripple stacks have been tested against an independent mesh (Supplementary Note 3 of the paper): the stack contact agrees to a few tenths of a per cent; the modelling choice of vertical offsets is not independently tested. The Level 3 record is in `level3/` with its input arrays. Redact local file paths from the two session transcripts before the public release: `python3 tools/redact_transcripts.py` is a dry run that reports 4 replacements in each, and `--write` applies them (keeping `.orig` copies, which must not be committed). Decide the freeze date and record it in `CHANGELOG.md`.

## At the freeze (v1.0.0-rc1)

Rebuild the checksums. Exam B is already designed, its predictions are filed and its truth is sealed (`code/exam_B/`); its truth file is held outside the repository and is added after scoring. If the law, the simulator or the exam specification changes after the human review, the predictions must be re-filed (and the truth re-simulated if the simulator changed), as `code/exam_B/criteria.md` explains. Re-run every number that appears in the paper from a clean clone and record the commit and the environment against it. Review `docs/` and rewrite the working documents as a single reproduction note.

## At submission (v1.0.0)

Choose the licence (a permissive one for the code and CC BY 4.0 for data and documents is the usual split; the University of Manchester's technology-transfer office should be consulted first because of the spin-out interests) and add `LICENSE`. Make the repository public. As I understand Zenodo's GitHub integration it archives releases of public repositories, so the DOI can only be minted once the repository is public; check the current Zenodo documentation before relying on this. Link the GitHub account to Zenodo, switch the repository on there, and create the GitHub release `v1.0.0`; Zenodo then issues a version DOI and a concept DOI. Put the concept DOI in `CITATION.cff`, `.zenodo.json` and the paper's Code availability statement, and commit that change before tagging so that the archive carries it. If reviewers need access before the repository can be public, either add them to the private repository by username or upload a zip of the release by hand to a draft record on Zenodo, which reserves a DOI before publication; confirm the current options with Zenodo.

## Decisions needed from the authors

The repository name and owner on GitHub. The licence and who holds copyright. The author list in `CITATION.cff`. Whether a public tagged version is needed earlier than the science paper: the method paper on the Co-Scientist evaluation (target 27 January 2027) will refer to the exam, the scoring harness and the verification protocol, so it may need something citable before the science paper's deadline. The Co-Scientist scoring tools (`harness_scoring/` in the project's release candidate) are not in this repository because they depend on the Co-Scientist outputs and the terms for redistributing those are unresolved.
