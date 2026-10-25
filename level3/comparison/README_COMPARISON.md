# Operator's guide to the reimplementation comparison (do not give this folder to the reimplementer or the reviewer)

## What is in this folder

pairs.npz and configs.npz hold the reference answers for the deterministic tests. reference_stats.json holds the reference means and standard errors. reference_axes holds the reference per-axis winding results for thin discs (L = 16, 24, 32, 64 and 96). generate_reference.py regenerates the reference files from the reference simulator. compare_reimplementation.py does the comparison.

## Order of work

First run the Cursor reimplementation folder to completion in its own chat (see its INSTRUCTIONS.md). Only then run the review folder, in a new chat. Copy results.json, work/ and the exported chat transcript out of the reimplementation folder before doing anything else. Then run, from any directory: python3 compare_reimplementation.py path/to/results.json

## How to read the output

Section A and B are exact tests. A handful of tangent-contact mismatches (up to 0.2 per cent of pairs, one configuration per set) are tolerated; anything beyond that means the two contact tests or the winding rules genuinely differ, and the pairs that differ should be listed and examined first. Section C compares means: pass means agreement within three combined standard errors plus an allowance of 0.5 per cent (0.2 per cent for the thin-disc plateau references, whose own uncertainty is small). Section D checks that the reimplementer's numbers are internally consistent; a failure here is a reporting error in their results, not a physics disagreement. Section E is diagnostic only.

## Decision logic

If everything passes, the thin-disc plateau of about 0.977 is reproduced by a code written from the specification alone. That makes a coding error in the reference simulator an unlikely explanation of the 1.6 per cent offset from the published 0.9614; the remaining explanations are a bias in the published extrapolation or a defect shared by both codes because the specification contains it (for example a convention in the specification that differs from the published model). The per-axis runs already exclude the choice of crossing criterion as the cause, since single-axis and any-direction criteria converge on the same limit. Report the offset as a discrepancy with the published value, say that two independent codes agree on it, and do not describe the code as validated against the literature.

If the thin-disc statistics fail but pairs and configurations pass, the contact tests agree and the difference lies in insertion, neighbour search, periodic images or the wrapping criterion. Section E then tells which: if the reimplementer's any-direction value is lower than the reference by about 1.5 per cent, towards 0.9614, the reference may have a wrapping or image defect; if the reimplementer's single-axis values match the reference single-axis values but the any-direction values differ, the definition of the criterion differs. Resolve by reading both codes at the insertion and wrapping routines before changing anything.

If a pair or configuration test fails, find the first differing item, check it by brute force, and decide which code is wrong. Do not average, tune or choose by which result looks closer to the literature.

If the cylinder or polydisperse cases fail while the thin-disc cases pass, the contact test for thick bodies, the log-normal draw or the aligned-domain generator differs; compare the generators' measured moments first (the reimplementer was asked to record P2 for kappa = 6.312 in NOTES.md).

Whatever the outcome, record it in 10_Release_Manifest.md and the provenance ledger, with the transcript and the reimplementer's code checksummed.


## Additions of 25 Oct 2026 (zero-thickness squares, triangles and oblate spheroids of aspect ratio 1000)

The reimplementation package now also asks for these three particle types (SPEC.md sections 1, 2, 5 and 7). Reference answers are in shapes_pairs_ref.npz and shapes_configs_ref.npz, produced by generate_reference_shapes.py (it needs literature_checks/sim2e.py and poly_test.py); statistical references are the entries square_*, triangle_* and spheroid1000_* in reference_stats.json, taken from the runs in literature_checks/shape_fss.jsonl and ell_fss.jsonl (different seeds from anything the reimplementer sees). compare_reimplementation.py prints them as sections A2, B2 and C2 after the original sections; the raw-array consistency check D now includes the plate cases.

How to read them. The reference plateau values are about 0.879 (squares), 0.750 (triangles) and 0.979 (spheroids of aspect ratio 1000), against the published thin-plate values of 0.8647 and 0.7295 for the first two and the 1.27 that follows from Garboczi et al. for the third. If the reimplementer reproduces the reference values, a coding error in the reference simulator becomes an unlikely cause of the 1.6, 1.7 and 2.9 per cent offsets from the published values, and of the difference from Garboczi (the polygon test is a different algorithm from the disc test, and the spheroid test is a different support function from the cylinder test). If the pair and configuration tests pass and the statistics fail, read the insertion and wrapping routines first, as for the discs. If the reimplementer lands on the published values for squares and triangles but not for discs, or the reverse, that is a finding about the specification or the contact test of that shape, and both codes must be read before anything is concluded. Do not average, tune or pick by which result looks closer to the literature.
