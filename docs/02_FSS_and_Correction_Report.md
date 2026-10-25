# Finite-size scaling and the 1.015 correction

Date: 24 Oct 2026. Scripts: verification/fss_thin.py, fss_thin_large.py, fss_spheres.py, vex_check.py, binary_check.py. Raw results: fss_thin.json, fss_L64.json, fss_L96.json, fss_L128.json, fss_spheres.json, binary_check.json.

## Question

The simulator's thresholds are divided by 1.015 before being compared with the law, and the law's thin-disc limit is anchored to the Yi and Tawerghi (2009) reduced critical density of 0.9614. The question was whether 1.015 is a finite-size effect that vanishes in large boxes, or a persistent offset.

## Method

Zero-thickness discs of unit radius are inserted one at a time into a periodic cube of side L, and the threshold is the number of discs present when the first cluster with non-zero winding appears, converted to the reduced density eta = (4 pi / 3) N / L cubed. Box sizes from 8 to 128 radii were run, with 300 realisations up to L = 48, 150 at L = 64, 60 at L = 96 and 30 at L = 128.

## Result

The mean reduced density is 1.0032 at L = 8, 0.9851 at 12, 0.9770 at 16, 0.9807 at 24, 0.9772 at 32, 0.9756 at 48, 0.9782 at 64, 0.9783 at 96 and 0.9761 at 128, with standard errors of 0.001 to 0.008. Beyond L = 32 there is no trend: the weighted mean of the five boxes from L = 32 upward is 0.9770 plus or minus 0.0006 (0.9769 without L = 32), which is 1.6 per cent above 0.9614. The earlier extrapolation that gave 0.9738 plus or minus 0.0020 assumed a power-law approach with exponent 1/nu = 1.141 and was wrong in assuming that a trend still existed; the larger boxes show a plateau. The offset is therefore not a finite-size effect in the ordinary sense. The empirical factor of 1.015 describes the plateau to within 0.1 per cent (the plateau ratio is 1.0162 plus or minus 0.0006), so using 1.015 does not introduce error beyond the 0.1 per cent level.

## Controls

The contact test is exact. A forty-million-pair Monte Carlo estimate of the mean excluded volume of two randomly oriented unit discs gives 9.8715 plus or minus 0.0024 against the analytic pi squared of 9.8696, so the offset cannot come from a missed or spurious disc-disc contact.

The same pipeline applied to spheres gives 0.3441 plus or minus 0.0009 at L = 32 against the accepted 0.3418 (volume fraction 0.28954, Lorenz and Ziff), a plateau 0.7 per cent high from L = 16 to L = 32. The wrapping criterion therefore adds a small positive offset even for spheres, but not the 1.6 per cent seen for discs.

A binary mixture check against the statement in Yi and Tawerghi that a 50:50 mixture of two disc sizes percolates within 1 per cent of the equal-disc value, with eta defined through the mean cubed radius, gave +1.4 per cent (standard error about 0.4) for a radius ratio of 0.5 and +0.9 per cent (about 0.4) for 0.25, in L = 40 boxes with 80 realisations. This is marginally above their statement, within about two to three standard errors, and is not a clean confirmation.

## Interpretation

## Per-axis winding and the estimator hypothesis (added 24 Oct, evening)

The first explanation offered in the earlier draft was that the any-direction winding criterion is a different estimator from the spanning probability used by Yi and Tawerghi. That was tested directly with verification/axes_check.py, which records for every realisation the first insertion count at which a closed path winds along each axis separately (single-axis criterion) as well as in any direction. Thin discs, box sides 16, 24, 32, 64 and 96, with 400, 300, 300, 100 and 40 realisations. The mean over the three axes of the single-axis reduced density is 1.0306, 1.0100, 1.0023, 0.9867 and 0.9827 respectively, falling monotonically; the any-direction value is flat at 0.976 to 0.981 over the same range, as before. A fit of the single-axis values to eta_inf + a L^(-1/nu), with nu fixed at the three-dimensional percolation value 0.8765, gives eta_inf = 0.9756 plus or minus 0.0012 with chi-squared 1.2 for three degrees of freedom; leaving the exponent free gives 0.973 plus or minus 0.005 and an exponent of 1.0 plus or minus 0.2, consistent with 1/nu = 1.14 (verification/axes_fss_fit.py and its log).

The two estimators therefore converge on the same limit, 0.9756 for the single-axis criterion and 0.9776 plus or minus 0.0009 for the any-direction plateau, which is what one expects: the choice of crossing criterion changes the finite-size shift but not the thermodynamic limit. This weakens the first explanation considerably. A difference in estimator (mean first-wrap count against fifty per cent crossing probability, single axis against any direction) cannot by itself produce a persistent 1.5 per cent offset in the limit, because all sharp criteria share the same limit. What remains of the estimator hypothesis is only a possible finite-size effect in the published extrapolation, which depends on how Yi and Tawerghi extrapolated and cannot be checked from their summary.

## Interpretation

Two explanations remain open, and a third has become less likely. The published value may carry its own bias, since the same paper notes that Garboczi and co-workers' extrapolation for the oblate limit is considerably higher, and its extrapolation procedure is not reproduced here. The simulator may have a subtle defect that the tests do not cover, for instance in the geometry of near-tangent contacts or in the cell-list handling at the box sizes used; the sphere control (+0.7 per cent at L = 32, not yet repeated with the single-axis fit) argues against a gross defect in the cluster and winding code but not against a small one. The estimator explanation is disfavoured by the per-axis result above. The independent reimplementation is the decisive check: if it obtains 0.976 to 0.977, the offset is a property of the stated geometry and the published value needs explanation; if it obtains 0.9614, the author's code has a defect that must be found.

## Consequences for the science paper

The law is calibrated to simulation divided by 1.015, which equals the Yi and Tawerghi scale. The paper should state this, give the raw plateau (0.977) alongside the benchmark (0.9614), and report threshold ratios, which are insensitive to the factor, wherever possible. If the independent reimplementation obtains the same plateau, the offset is a property of the stated model as simulated (not of the crossing criterion) and should be reported as a discrepancy with the published value; if it obtains 0.9614, the author's code has a defect that must be found.

## Status of test T8

T8 (extrapolated eta_infinity equal to the benchmark within error) fails, and the failure is a real, reproducible and reproducible offset whose cause is unresolved (estimator excluded as a sole explanation, see above) of 1.6 per cent. It is recorded as a known discrepancy and is not to be deleted from the suite.


## Note added 25 Oct 2026

The thin-plate offset discussed above was extended on 25 Oct to squares, triangles and spheroids and to the published open-face procedure; see 15_Shape_Ripple_FSS_Checks_25Oct.md. Where this report and that document differ, the later document is the author's current position.
