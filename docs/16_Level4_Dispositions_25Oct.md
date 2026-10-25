# Level 4 (Grok, Run 2): dispositions after independent reproduction

Date: 25 October 2026. Reviewer: Grok 4.7 in Cursor (findings in `04_Simulation_and_analysis/verification_runs/Level4_Grok_run2/findings/FINDINGS.md`). Reproduction by Claude: every finding was re-run from Grok's own scripts (`test_law_docs.py`, `test_geometry.py`, `test_followup.py`) against a copy of the repository `code/` that was byte-identical to the review package, before any finding was accepted. Outputs are kept beside this note as `L4_repro_law_docs.out`, `L4_repro_geometry.out` and `L4_repro_followup.out`.

## Summary

All five findings reproduce exactly as described. None changes a reported threshold, a plateau, an exam error or a pre-registered test result. One (F2) corrects a description that is wrong in the ledger. Two (F3, F5) called for small corrections, to code and to wording; both have been made. One (F1) is an accepted, already disclosed limitation with one wording slip. One (F4) is an accepted limitation inside the stated box condition, and it retires a worry in the self-audit.

## F1. Thin-plate offset and circularity: confirmed, accepted limitation

The weighted disc mean for L of 32, 48, 64, 96 and 128 reproduces as 0.9770 plus or minus 0.0006, and 0.9769 for L of 48 and above. The face-convention over-marking reproduces as 7 of 1172 (0.60 per cent) and acts in the wrong direction to explain the offset. The winding number agrees with a 3x3x3 supercell on all four positive draws and is invariant to translation, a full-box shift and a 90 degree rotation. The offset is therefore in the data and is not produced by any of the tested detector defects. The circularity is real: c0 is 1.5 times 0.9614 by construction and every fitted threshold is divided by 1.015, so agreement of the law with the published thin-disc value is not evidence and must not be cited as confirmation. Actions: state this sentence in the main paper; correct document 02, where "four largest boxes" should read "five boxes from L = 32 upward"; the offset itself remains an open question for the human assessor and for the Level 3 convention check.

## F2. Kell(a) is not the arc-length ratio: confirmed, documentation error

The table in `test_law_docs.py` reproduces exactly. The formula in the code and the main text is the inverse excluded-volume ratio of an elliptical footprint with semi-axes k and 1, not the arc-length ratio of a sinusoid. The code is not wrong and the exam does not use the ripple factor. Actions: correct the ledger line; in the paper describe Kell as an elliptical-footprint ansatz fixed by geometry apart from the fitted damping CRIP; report the ripple fit quality (chi-squared 241 on 29) next to it, as it is the weakest fit in the law.

## F3. Polygon contact: confirmed, corrected in v0.1.1, no effect on results

Exactly parallel polygons returned no contact, and a vertex or edge lying in the other plane was counted only when the rest of the polygon lay on the positive side of the plane normal. Both reproduce (coplanar overlap code False against independent True; one-vertex case True on the positive side and False on the negative side). Random non-parallel pairs agreed with the independent chord test on 2000 of 2000 draws (703 contacts), so the square and triangle plateaus stand. Accepted by the corresponding author on 25 October 2026 and corrected the same day: the sign test in `poly_interval` is symmetric, the parallel case is decided by a separating-axis test in the common plane (`poly_coplanar_overlap`), and `run_once` refuses a box smaller than four times the largest bounding radius (F4). Evidence that the correction changes nothing else: 400,000 random square pairs, 0 differences between old and new; threshold runs for squares, triangles and discs (six realisations each, box of 16 radii) identical; the new `tests/test_polygon_contact.py` (4,000 random coplanar pairs per shape against a polygon-clipping oracle, vertex and edge cases from both sides and with roles swapped, 3,000 random pairs per shape against an independent chord test) passes; the quick regression suite (14 of 14), the larger-box test and the sim2/sim2e equivalence (36 configurations) pass. The first, partial patch (`level4/L4_F3_proposed_patch.diff`) is kept as a record. A vertex touching the other polygon with the rest on one side now counts as contact from both sides, a measure-zero event for which the reviewer's independent test says no contact; only symmetry was required.

## F4. Single minimum image: confirmed, accepted limitation inside the stated condition

The two-large-disc case reproduces (simulator returns no percolation, supercell returns a cluster of 2). The stated condition, box side greater than four times the largest bounding radius, is the right one and every reported campaign meets it. For exam cases F5 and F6 at L = 32 with radii capped at 6, the margin L/2 minus twice the largest bounding radius reproduces as 4.15 and 4.14, so no second-image contact is possible. The self-audit phrase "met only narrowly" for these two cases is retired. The condition is now stated in the paper's Methods and enforced in the code (`run_once` raises an error for a violating box).

## F5. Sphere-control sentence: confirmed, wording to fix

The plateau for L of 20 and above reproduces as 0.3422 plus or minus 0.0013 and the six-point 1/L fit as 0.3386. The exact-sphere campaign gives 0.3437, 0.3443 and 0.3441 plus or minus 0.0009 at L of 16, 24 and 32. The sentence should report the plateau and the exact-sphere value next to the extrapolation, and should claim only that the winding code reproduces the sphere threshold at the 0.3 to 0.7 per cent level, not a low bias of one per cent.

## What remains for human judgement

Whether the offset of 1.6 to 2.9 per cent is a shared specification convention (Level 3 shows it is not a coding error in `sim2e.py`); the unverified items the reviewer lists (the methods of Yi and Tawerghi, field-by-field ripple scatter, prolate resimulation, the experimental transcription); and the independent human review of the simulator.
