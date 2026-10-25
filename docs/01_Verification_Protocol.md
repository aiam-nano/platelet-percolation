# Verification protocol for the percolation simulator and the calibrated law

Status: DRAFT for approval (Claude drafted, 24 Oct 2026). It becomes binding once the operator approves it and the freeze date below is recorded in PROJECT_STATUS.md.

## 1. Purpose and scope

The science paper rests on three things that must each be shown to be correct before any number from them is published: the simulator that produces the "truth" values (sim2.py and its helpers), the calibrated law (finallaw2.py and mylaw.py), and the way the law was scored against held-out cases (the sealed twelve-case exam). Every one of these was written and run by the same AI system (Claude), which is also the analyst. That is a conflict of interest and a single-point-of-failure, so this protocol assumes the code is wrong until independent routes agree, and it records which checks are independent of the author and which are not.

## 2. Levels of evidence, from weakest to strongest

Level 0 is "the code runs and gives plausible numbers". It is not evidence and is never cited as such.

Level 1 is the regression suite written by the author (verification/test_suite.py). It checks the code against analytic or brute-force answers the author could compute independently of the implementation under test: exact disc contact against a brute-force point-sampling test, GJK cylinder contact against the analytic disc in the thin limit, mean excluded volume against the Isihara formula for cylinders and against pi squared for discs, completeness of the cell list, the wrapping-cluster detector against a brute-force union-find, the orientation generators, and the thin-disc finite-size trend. The suite passes 14 of 15 tests; the failing one (T8) is explained in 02_FSS_and_Correction_Report.md. Level 1 catches implementation slips but shares the author's understanding of the problem, so it cannot catch a misunderstanding.

Level 2 is comparison with other groups' published numbers: the thin-disc threshold of Yi and Tawerghi (2009), their binary-mixture statement, the sphere threshold of Lorenz and Ziff, and the Garboczi et al. oblate-ellipsoid limit. Published values can themselves disagree (Yi and Tawerghi report that their disc value is noticeably lower than the Garboczi extrapolation), so Level 2 gives a plausibility window, not a proof.

Level 3 is independent reimplementation: a second AI system, or a person, writes the simulator from the specification alone (04_Reimplementation_Spec.md) without seeing sim2.py, and both codes are run on the shared cases in the comparison harness. Agreement within the stated tolerances is the main evidence the simulator is right. This level has not yet been carried out and needs a second-system account.

Level 4 is adversarial review: a different AI system, and ideally a colleague, reads the code and the self-audit (05_Self_Audit_and_Review_Packet.md) and tries to break it. Not yet carried out.

## 3. What must pass before the science paper is drafted

Level 1 must pass in full or each failure must carry a written, accepted explanation. Level 2 checks must lie within the published uncertainties or the discrepancy must be stated in the paper. Level 3 must show that the independent code reproduces the benchmark cases to within the sampling error plus the finite-size tolerance in the spec, and the three non-trivial geometries (finite-thickness cylinders, aligned domains, rippled sheets) to within 3 per cent. Level 4 findings must be closed or listed as limitations. A claim in the paper that rests only on Level 1 must say so.

## 4. Freeze rule

The simulator, the law and the exam truth values are frozen when the operator records a freeze date. After that date any change to sim2.py, align.py, finallaw2.py, mylaw.py or the exam cases is a new version with a new tag, and every figure in the paper must name the version that produced it. The exam truth file is checked against its SHA-256 (truth_sealed.sha256) every time the exam is scored. The law's constants were fitted before the exam was run; they are not to be refitted to the exam. If a constant is refitted after seeing exam results, the exam is demoted to a calibration set and a new sealed set must be generated.

## 5. Uncertainty reporting rule

Every simulated threshold is reported with its standard error, the system size, the number of realisations and the correction factor applied. The finite-size factor 1.015 is a benchmark-tied empirical factor, not a theorem, and is quoted in every table that uses it. Percentage errors between law and simulation are quoted against the combined sampling error, and any agreement tighter than the sampling error is not claimed.

## 6. Handling of post hoc choices

Three choices made after seeing exam truth are disclosed and not to be described as predictions: the mapping from the prompt's director order parameter to the law's operational local order parameter (S_local = 0.665 S_dir squared, measured on the generator), the choice of domain size d = 6 radii instead of the naive d = 3, and the conversion of the law's nominal volume fraction to the exam's occupied volume fraction through 1 minus exp(minus phi). Each must appear in the Methods with its justification. A fresh sealed set generated after the freeze would remove the need for this disclosure and is recommended if time allows.

## 7. Roles

The author (Claude) writes code, runs Level 1 and Level 2, and prepares the packets. The operator approves this protocol, supplies a second system for Levels 3 and 4, and arranges a human reviewer. A second assessor (not the author) scores any hypothesis-quality battery for the method paper.

## 8. Record

Each verification run is logged with date, git or checksum identifier of the code, the seed and the output file. The checksums are listed in 10_Release_Manifest.md. Failures are never deleted from the record.


## Approval and implementation record (24 Oct 2026)

Approved by the principal investigator on 24 Oct 2026 as written. Levels 1 and 2 are run and recorded by Claude. Levels 3 and 4 are run by Grok 4.7 in Cursor, as two separate co-work folders: 04_Simulation_and_analysis/reimplementation_cursor (Level 3; answers and comparison script held back in reimplementation_HELDBACK) and 04_Simulation_and_analysis/review_cursor (Level 4, run afterwards in a new chat). The science paper is geometry-only (no tunnelling or electronic properties). The alignment prefactor is refitted to 0.588 (finallaw3.py) and a fresh sealed exam is generated after the freeze.
