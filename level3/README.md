# Level 3: independent reimplementation (Run 1, 25 October 2026)

A second AI system (Grok 4.7 in the Cursor editor) rebuilt the simulator from `package/SPEC.md`, which contains no source code and no reference result, in an empty workspace. Its outputs were compared with the original code by `comparison/compare_reimplementation.py`, which it never saw. The criteria were fixed before the run: at most 0.2 per cent of pairs and at most 1 of 20 configurations may differ, and means must agree within three combined standard errors plus 0.5 per cent (0.2 per cent for the thin-disc plateau references).

`comparison_output.txt` is the comparison script's output: OVERALL PASS. No mismatch in 5 x 4,000 contact decisions or in 6 x 20 configurations; all 14 statistical comparisons pass; the largest difference is 0.82 per cent.

Contents. `package/` is what the reimplementing system was given; its input arrays are in `package/inputs/` (4.3 MB), and the same arrays together with the reference answers are in `comparison/` (`pairs.npz`, `configs.npz`, `shapes_pairs_ref.npz`, `shapes_configs_ref.npz`). `tests/verify_level3_inputs.py` regenerates every one of them with `comparison/generate_reference.py` and `generate_reference_shapes.py` and checks that they are identical bit for bit (the statistical references in `reference_stats.json` are long simulations and are not rerun). `comparison/` holds the comparison script, the reference generators and the reference statistics, which were withheld from the reimplementing system. `work/` is the reimplementing system's own code (`core.py`), its self-tests, notes and report. `results.json` is its returned result file and `cursor_instructions_execution_process.md` is the transcript of the session.

The session transcript contains local file paths from the operator's computer; remove or redact them before the public release.

Limits. The specification was written by the same author as the original code, so a convention shared by both would not be detected. Ripple was not part of the package.
