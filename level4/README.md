# Level 4: adversarial review (Run 2, 25 October 2026)

Grok 4.7 in the Cursor editor reviewed the code and the audit packet in a separate workspace (no reimplementation, no held-back references). `FINDINGS.md` is its report: five findings, F1 to F5. `review_work/` holds its three reproducing scripts and notes. They were written to run from a folder that contains a copy of this repository's `code/` directory named `code_copy`, as it was at version 0.1.0; to run them now, copy `code/` to `review_work/code_copy` (the polygon behaviour in `test_geometry.py` will then differ from the transcript, because version 0.1.1 corrects it).

The findings were not accepted on the reviewer's word. Each reproducing script was run again by the AI agent that wrote the code, against an identical copy of the version 0.1.0 code, and all five reproduced; the outputs are `review_work/L4_repro_*.out`. The dispositions are in `../docs/16_Level4_Dispositions_25Oct.md`. The one code defect (F3, polygon contact) was corrected in version 0.1.1; `L4_F3_proposed_patch.diff` is the first, partial patch (symmetric sign test only). The final change also handles exactly parallel plates and is covered by `tests/test_polygon_contact.py`.

The session transcript contains local file paths from the operator's computer; remove or redact them before the public release.
