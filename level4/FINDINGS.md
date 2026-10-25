# Findings

The disc, cylinder and spheroid contact tests, the cell list in the boxes this code actually uses, and the winding detector all behaved correctly on the cases I could check by running them. The thin-plate plateaus this simulator reports, about 0.977 for discs, 0.879 for squares and 0.751 for triangles, are what the packaged runs contain, and I did not find a defect that moves them down to the published 0.9614, 0.8647 and 0.7295. The law's thin-disc limit is not an independent confirmation of 0.9614: the constant is that number, and every fitted threshold was divided by 1.015 first. The ripple factor is not the arc-length ratio the ledger says it is. The exam arithmetic and the "22 of 24, with two named failures" tally are correct. Polygon contact is wrong for parallel plates and for one class of vertex and edge contacts; that class has no effect on the random-orientation plateaus.

## F1. The thin-plate offset is real in these data, and the law's agreement with 0.9614 is built in

**Where.** `code/sim2.py` line 325 (the threshold is the count of the inserting particle). `code/mylaw.py` line 3 (`c0 = 1.4421`), used by `code/finallaw3.py` line 15. `code/fitlaw.py` line 5 and `code/calibrate_all.py` line 22 (`FS = 1.015`). The numbers are in `docs/02_FSS_and_Correction_Report.md` line 15 and `docs/15_Shape_Ripple_FSS_Checks_25Oct.md` line 15.

**What is wrong or doubtful.** Two separate things.

The plateaus sit above Yi and Tawerghi, Phys. Rev. E 79, 041134 (2009). From the packaged outputs, the 24 October disc campaign gives a weighted mean of 0.9770 ± 0.0006 over L = 32, 48, 64, 96 and 128, which is 1.0162 times 0.9614. The later shape run, restricted to L ≥ 32, gives 0.9767 ± 0.0008 (discs, factor 1.0159), 0.8794 ± 0.0007 (squares, 1.0170) and 0.7505 ± 0.0006 (triangles, 1.0287). Square/disc is 0.9004 against a published 0.8994. Triangle/disc is 0.7683 against 0.7588, about 1.2% high, so the offset is not one common factor. I did not find a bug that produces it. Disc excluded volume from the contact test is 9.879 ± 0.019 against π² = 9.870 (quick suite, T4). The cell list misses no pair in direct comparisons at cell counts 8, 12 and 15. The winding number matches a 3×3×3 supercell on four percolating draws and on one hand-built chain. The GJK sphere plateau of the spheroid code, L ≥ 20, is 0.3422 ± 0.0013 against 0.3418, so the detector is not globally 1.6% high. An open-boundary reproduction with hard faces, extrapolated in 1/L, still lands at 0.9745 (median) and 0.9728 (mean) for discs. The face test marks a disc that crosses the infinite plane of a face even if the crossing lies outside the face square; on 4000 random discs that happens for 7 of 1172 face marks, 0.6%, and that error would lower a threshold, not raise it. I could not read the methods section of Yi and Tawerghi (it is paywalled), so I cannot say their extrapolation matches `span_open.py` beyond the abstract.

The second point is circularity, and it does affect a claim. `c0 = 1.4421` is 1.5 × 0.9614, so the thin-disc limit of the law is 0.9614 by arithmetic. Every simulated threshold passed to the fit is divided by 1.015, which is the disc offset and not the square or triangle offset. Absolute agreement of the law with the published thin-disc number is therefore not evidence. Ratio predictions, the gap tests, and the exam cases do not use that agreement; they compare the law with this simulator. The author says all of this in the ledger and in document 02. The sentence in document 02 that says "four largest boxes" is the only slip: the quoted 0.9770 ± 0.0006 is the weighted mean of the five boxes from L = 32 upward. Dropping L = 32 changes it to 0.9769.

**Reproducing case.** `review_work/test_followup.py`, from `review_work`:

```
python test_followup.py
```

The large-box block prints L = 32: 0.9772, 48: 0.9756, 64: 0.9782, 96: 0.9783, 128: 0.9761, and the weighted mean of those five is 0.9770 ± 0.0006. The face block prints 7 chord misses out of 1172 extent marks. The winding block prints four pairs of equal thresholds (80, 91, 99, 90 on the seed used) and zero mismatches. `review_work/test_law_docs.py` prints the thin-disc limit as 0.961397. `review_work/test_geometry.py` prints zero missed neighbour pairs. The quick suite T4 prints the excluded-volume comparison.

**Effect.** About +1.6% on the thin-disc threshold, +1.7% on squares, +2.9% on triangles, relative to the 2009 values. On the law, replacing 1.015 by 1.0159 moves every absolute threshold by 0.09%. The circularity removes the thin-disc literature number as an independent check of the absolute scale. It does not move the exam percentages, which are ratios against this simulator.

**Confidence.** High that the offset is in the data and that the tests above did not produce it. Medium that no such bug exists: an independent implementation of the same estimator would decide it, and this review was not allowed to look at one. High that the circularity is as described.

**Disposition.** Needs human judgement for the offset. The circularity is an accepted limitation, already disclosed, and it should not be cited as a confirmation of 0.9614.

## F2. Kell(a) is not the arc-length ratio

**Where.** `code/finallaw3.py` lines 6–9, and `code/finallaw2.py` lines 5–8. The claim is `docs/03_Provenance_Ledger.md` line 15: "The elliptic-integral factor Kell(a), the arc-length ratio of a sinusoidal ripple of slope parameter a, is a closed-form geometric quantity and has no free constant."

**What is wrong.** The arc-length ratio of one period of z = a sin θ, relative to the flat length, is

α(a) = (2/π) ∫₀^{π/2} √(1 + a² cos² θ) dθ = √(1+a²) · (2/π) E(a²/(1+a²)),

where E is the complete elliptic integral of the second kind. The code estimates α by a uniform average of the integrand and then sets kx = 1/α and

Kell(a) = π / (2 · kx · E(1 − kx²)) = α π / (2 E(1 − 1/α²)).

That is not α. At a = 1.2, α = 1.2965 and Kell = 1.4578, a ratio of 1.124. The excess Kell − 1 is 1.51 times α − 1 at a = 0.5 and 1.60 times α − 1 at a = 3. The same expression is in both law files, and CRIP was fitted against it, so the in-sample ripple simulations and exam cases F10 and F11 still match this formula. What fails is the statement that the factor is the arc-length ratio and has no fitted content in its a dependence. The ripple fit itself is poor on its own terms: χ² = 241 on 29 degrees of freedom in `verification/calibrate_all.log`.

**Reproducing case.** `review_work/test_law_docs.py`:

```
python test_law_docs.py
```

The Kell block prints, among other rows,

```
a=1.2  arc_quad=1.296484  arc_ellip=1.296484  Kell=1.457783  Kell/arc=1.124412  k=0.771317
```

**Effect.** Keeping CRIP = 0.647 and replacing Kell by α at a = 1.2, λ = 1 changes the ripple ratio from 1.2106 to 1.136, which is 6.5% relative. A refit of the one constant CRIP would absorb most of that, because (Kell − 1)/(α − 1) is nearly constant; the leftover shape error is about one percent of the exam ratio. The exam scores, which used this Kell, do not change. Extrapolation in a, and any claim that the a dependence is fixed by geometry, do change.

**Confidence.** High. The identity is a few lines and the numerical comparison is exact.

**Disposition.** Fix the formula, or stop calling it the arc-length ratio and treat its a dependence as part of the fit.

## F3. Polygon contact misses parallel overlaps and one-sided vertices and edges

**Where.** `code/literature_checks/sim2e.py`, `poly_contact` lines 275–295 and `poly_interval` lines 256–268. The parallel branch is the `return False` at lines 285–286. The sign test is line 258, `(s[k] <= 0) != (s[k2] <= 0)`.

**What is wrong.** Two defects.

Exactly parallel polygons always return no contact, including two coplanar squares that overlap. The disc code treats that case; the polygon code does not. Random orientations hit exact parallelism with probability zero, so the square and triangle plateaus are unaffected. A run with aligned plates would never percolate.

A vertex or an edge lying in the other plane is counted only when the rest of the polygon is on the positive side of the plane normal. I placed one vertex of a square at the origin, inside the other square, with the other three vertices at z = +1.77. The code returns contact. Flipping the normal, so the same vertex is still at the origin and the body lies at z = −1.77, the code returns no contact. An edge through the other square does the same: contact on the positive side, no contact on the negative side. The test is the comparison `s <= 0`, which treats a point in the plane as "already negative", so a crossing is seen only when the other endpoint is positive.

Random non-parallel pairs are fine: 2000 random squares gave 0 disagreements with an independent chord test (703 contacts).

**Reproducing case.** `review_work/test_geometry.py`:

```
python test_geometry.py
```

The polygon block prints `coplanar overlap code False independent True`, the vertex coordinates with one vertex at the origin, `vertex-in-plane, polygon on + side: code True`, and `vertex-in-plane, polygon on - side: code False`. The edge block prints `edge-in-plane code True` and `edge other side ... code False`.

**Effect.** None on the reported random-orientation thresholds, at the precision of those runs. A perfectly aligned polygon simulation would be entirely wrong.

**Confidence.** High. The negative-side miss follows from the sign test and from the coordinates.

**Disposition.** Fix. The random-orientation plateaus can stand.

## F4. One minimum image is recorded, so a pair that overlaps twice is not a wrap

**Where.** `code/sim2.py` lines 304–307. The same pattern is in `sim2e.py` at the contact call inside `first_perc`.

**What is wrong.** Contact is tested on the single minimum-image vector. If two particles also overlap through a second image, that second contact is never inserted, so the union-find never sees the offset that would wind. Two coplanar discs of radius 4 in a box of side 10, centres 3 apart, overlap on the near image (distance 3) and on the far image (distance 7). The simulator returns no percolation. A 3×3×3 supercell returns a cluster of 2.

The author states the condition L greater than four times the largest bounding radius. It is the right condition, and the reported campaigns meet it. For the exam cases that the self-audit flags as close, F5 and F6 at L = 32 with radii capped at 6, the largest bounding radius in a draw of 200,000 particles is 5.93. Half the box is 16 and twice that radius is 11.85, so two particles at the cap still miss by 4.1. No second-image contact is possible there. The phrase "met only narrowly" overstates the risk for those two cases.

**Reproducing case.** `review_work/test_geometry.py` prints `two large discs sim -1 supercell 2`. `review_work/test_followup.py` prints the F5 and F6 margins, `L/2-2max` about 4.15.

**Effect.** None on the thin-disc plateau, the shape plateaus, or exam cases F1–F9 and F12, because those boxes satisfy the condition. The failure is real for a box smaller than twice the sum of the two radii.

**Confidence.** High.

**Disposition.** Accepted limitation inside the stated condition. The F5/F6 worry in the self-audit can be retired for these caps.

## F5. The sphere-control sentence quotes an extrapolation the plateau does not support

**Where.** `docs/05_Self_Audit_and_Review_Packet.md` line 27, and `code/literature_checks/sphere_fss.py` lines 9–11. The fit is in `code/literature_checks/sphere_fss.log`.

**What is wrong.** The sentence says a spheroid run at aspect ratio 1, extrapolated against 1/L, gives 0.3386 ± 0.0019 against the accepted 0.3418. The fit is what the log contains, and its χ² is 4.8 on 4 degrees of freedom, but the points are not a clean 1/L line: L = 40 is 0.3397 and L = 56 is 0.3417. The weighted mean of L ≥ 20 is 0.3422 ± 0.0013. The separate exact-sphere campaign at L = 32 is 0.3441 ± 0.0009, which is the +0.7% document 02 reports, and it is flat from L = 16 (0.3437) to L = 32 (0.3441). Quoting 0.3386 as the control makes the detector look about 1% low. The plateau of the same runs is consistent with 0.3418.

**Reproducing case.** `review_work/test_law_docs.py` prints `weighted mean L>=20: 0.3422 +- 0.0013` and `fit all six: 0.3386`. `review_work/test_followup.py` prints the exact-sphere means, including L = 32: 0.3441 ± 0.0009.

**Effect.** No change to a threshold. It changes how much the sphere control can be asked to carry: it supports the winding code at the 0.3% level of that plateau, and it does not show a 1% low bias.

**Confidence.** High on the arithmetic.

**Disposition.** Fix the sentence so it reports the plateau next to the extrapolation.

## Checked and fine

- Quick regression, `python3 test_suite.py quick` in the verification copy: 14 of 14 passed in about 80 seconds. T8 is not in the quick set. The author says the full set fails T8; that is the offset in F1, and I did not rerun the long Monte Carlo.
- Disc contact. Parallel coplanar overlap is contact; an exact edge touch is not, because the test is strict `<`; a gap of 0.001 in radius is contact. Antiparallel normals behave as parallel. Perpendicular discs that cross are contact. A tilt of 10⁻⁶ uses the chord test; a tilt of 10⁻¹² with a plane gap of 10⁻⁸ is rejected by the parallel test. Particle swap: 0 disagreements in 2000 pairs (`test_geometry.py`).
- Cylinder support. Face-on, edge-on, a diagonal, and a direction 10⁻⁸ off the axis all matched a surface sample at half-heights 0.05, 0.5 and 2, with zero disagreements. Parallel cylinders at half-heights 0.001, 0.2 and 2: contact at gaps −0.05, −0.01, −10⁻⁴ and 0, no contact at +10⁻⁴, +0.01 and +0.05. A distance minimiser agreed on 185 random pairs, with 15 optimiser failures and 0 swap mismatches.
- Spheroid support is the analytic ellipsoid support. It matched surface sampling at half-heights 1, 0.1 and 0.001. At half-height 0.001, 0 of 1500 pairs disagreed with the disc test. Face-on, centre separation 0.0015 (overlapping slabs of half-thickness 0.001) is contact and separation 0.01 is not. Edge-on, separation 1.5 is contact and 2.1 is not.
- GJK against exact spheres, equal radii: 0 disagreements in 3000 pairs, including pairs closer than 10⁻⁶. The 60-iteration cap did not misclassify any pair I tried.
- Cell list. Against all pairs: 0 missed and 0 extra, for a monodisperse box (15 cells), a broad radius distribution (8 cells) and particles on the boundary (12 cells). The official T5 is vacuous whenever the cell is larger than twice the largest bounding radius, because a miss is then impossible; the author says so. The direct comparison covers the geometries T5 cannot fail.
- Union-find. `find` does not compress paths; it sums offsets along the parent chain. A chain of four discs that closes through two different images gives threshold 4 from the simulator and from the supercell. Shifting a percolating configuration by (1.3, −2.7, 4.4), by a full box, or rotating it by 90 degrees about z, leaves the threshold unchanged (`test_followup.py`).
- Rotations. Gaussian quaternions are the Haar measure on rotations, and the disc and cylinder only use the normal. The Fisher inverse CDF in `align.py` is the standard one. T7 passed, including κ = 4.733 giving S ≈ 0.50 and κ = 6.312 giving S ≈ 0.60, and the periodic Voronoi assignment.
- Polydispersity matrix. For a log-normal of number-mean radius 1, the eigenvalue ratio in `mylaw.py` equals the Gauss–Hermite integral of the Isihara kernel at full thickness T = 2/A, to four digits, at σ = 0.2, 0.3, 0.4, 0.5 and 0.8 and A = 8, 10, 20, 30, 60, 80 and 10⁶ (`test_law_docs.py`). The factors of 8, 4 and 2 are consistent once T = 2/A is used. This is not a bug.
- Ae = A/N. `phi(A=60, N=4)` equals `phi(A=15)` exactly, which is the right equivalent of the exam's single cylinder of thickness 4 × (2/60).
- Log-normal cap. For σ ≤ 0.5 and a cap at 6, the probability of a draw above the cap is about 7 × 10⁻⁵ and the second-moment bias is 0.2%. For σ = 0.8 and a cap at 5 the bias is large (the second moment falls by about 20%), and both the simulation and the registered prediction for that case renormalise the realised sample, so the comparison is internal.
- Exam scoring. With the documented shell values S = 0.1662 and 0.2394 and d = 6, `finallaw2` reproduces document 06 to the printed precision: errors +0.33, +0.41, −0.41, −0.05, +1.74, +1.60, +4.21, +1.90, +0.77 percent, mean absolute error 1.271%. `finallaw3` moves the three aligned cases to +1.07, +1.01 and +3.80, mean absolute error 1.085%. The ledger's "1.1, 1.0 and 3.8" matches. The shell values were not remeasured; they reproduce the published scores. The aligned cases were mapped after the exam was scored, as the ledger says, so they are not a blind test.
- Gap tests. Re-scoring the packaged outputs gives the stated 22 of 24. B2 is +18.7% and H1 is +23.6%, both outside ±8%, and every other case is inside its window, including TB5 at −3.2% inside its registered ±15%. The predictions are computed from the law and from the eigenvalue routine; they are not equal to the simulation results, so the central values were not fitted to these runs. The windows are the author's.
- Five fitted constants, not three: B_, κ, the alignment prefactor, γ, and CRIP. The ledger lists them and their χ². The alignment exponent was fixed at 2 and then released; the free fit is the one in `finallaw3`.
- Stacks, monodisperse cylinders and the G1 polydisperse cases agree with the law at the level the gap logs state.
- Ripple mesh. The reference contact (triangle–triangle) does not call the production contact (edges against the implicit surface). Both place the sheet on z = A sin(kx) with the same arc-length reparameterisation, which the document discloses. On the fine mesh the curved ratios are 0.9957, 0.9946, 0.9935, 0.9964 and 0.9941, i.e. −0.43% to −0.65%, and the production side is the one that misses pairs. The document's "0.4 to 0.7 percent" and the sign are right. The quoted band 0.9941–0.9964 leaves out the 0.9935 row; the prose range covers it.
- Open-boundary medians from `span_open.jsonl`: disc 0.9745, square 0.8743, triangle 0.7420. The document's statement that the mean and the median agree to 0.0003 is true for the squares and the triangles. The disc mean and median differ by 0.0017.
- Experimental CSV. All 16 rows match `finallaw3.phi` on the mid ratio and on the law's volume percent. Five mid ratios lie between 0.7 and 1.3. The extremes are 0.14 and 186. The arithmetic claim is right. The transcription from the papers was not checked.
- Equal-area polygons. The square of side √π and the triangle of side √(4π/√3) both have area π.
- Document 15's Garboczi sentence is the consistent one. The Wikipedia compilation gives occupied fraction 0.001275 and nV = 0.001276 at aspect 1000, hence (4/3)π n R³ = 1.276. 0.979 is 23% below 1.27, and 1.27 is 30% above 0.979. The packet's "about 30% higher" matches the second of those.

## Unverified

- Whether `span_open.py` matches Yi and Tawerghi's procedure. I have the abstract and the author's description, not the methods section. The face-extent convention over-marks by about 0.6%, which is too small, and in the wrong direction, to explain the offset.
- The per-field correlation of ripple test R2 is −0.08 on 10 fields. The scored test only checks the mean ratio. I did not pursue the field-by-field scatter.
- Prolate thresholds were not resimulated. The support function used for them passed the sphere limit and the sampling checks.
- The five thin-limit polydispersity residuals (the Br column) have no separate log in this folder. The reading that they are residuals, not full ratios, rests on the loss function and on G1.
- S was not recomputed with `shell_P2`. The documented values reproduce the scores, which is all the arithmetic needs.
- Experimental numbers were not checked against the papers. Several rows of the CSV already say the units look uncertain.
- `Kal` is negative for S = 1 and d = 1. That point is outside the fitted domain, which starts at d = 2, and I did not treat it as a finding.
