# Route to v1.0.0

Working note for the authors. Nothing here has been done unless it is marked done.

## Where things stand

v0.1.0 (25 October 2026) is a private, pre-freeze snapshot of the simulator, the law, the calibration, the regression suite, the verification scripts and the reproduction note. It has no licence and no DOI.

## Before the freeze (v0.x)

Settle the Level 4 findings and the human review, and fix or record each as an accepted limitation. Decide how to state the 0.4 to 0.7 per cent ripple-contact bias and whether to refit the ripple constant at a finer mesh. Test layered ripple stacks, or state that they are out of scope. Add the Level 3 record under `level3/`: the code-free specification, the cases, the input files, the validation script, Grok's code and report, and the comparison output. The held-back reference arrays and the comparison script can be published once Level 3 has been accepted. Decide the freeze date and record it in `CHANGELOG.md`.

## At the freeze (v1.0.0-rc1)

Rebuild the checksums. Generate a fresh sealed exam: commit its SHA-256 only, keep the truth file out of the repository until the predictions have been scored, and then add it. Re-run every number that appears in the paper from a clean clone and record the commit and the environment against it. Review `docs/` and rewrite the working documents as a single reproduction note.

## At submission (v1.0.0)

Choose the licence (a permissive one for the code and CC BY 4.0 for data and documents is the usual split; the University of Manchester's technology-transfer office should be consulted first because of the spin-out interests) and add `LICENSE`. Make the repository public. As I understand Zenodo's GitHub integration it archives releases of public repositories, so the DOI can only be minted once the repository is public; check the current Zenodo documentation before relying on this. Link the GitHub account to Zenodo, switch the repository on there, and create the GitHub release `v1.0.0`; Zenodo then issues a version DOI and a concept DOI. Put the concept DOI in `CITATION.cff`, `.zenodo.json` and the paper's Code availability statement, and commit that change before tagging so that the archive carries it. If reviewers need access before the repository can be public, either add them to the private repository by username or upload a zip of the release by hand to a draft record on Zenodo, which reserves a DOI before publication; confirm the current options with Zenodo.

## Decisions needed from the authors

The repository name and owner on GitHub. The licence and who holds copyright. The author list in `CITATION.cff`. Whether a public tagged version is needed earlier than the science paper: the method paper on the Co-Scientist evaluation (target 27 January 2027) will refer to the exam, the scoring harness and the verification protocol, so it may need something citable before the science paper's deadline. The Co-Scientist scoring tools (`harness_scoring/` in the project's release candidate) are not in this repository because they depend on the Co-Scientist outputs and the terms for redistributing those are unresolved.
