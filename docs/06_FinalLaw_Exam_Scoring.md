# Scoring the calibrated law (finallaw2) on the sealed twelve-case exam

Date: 24 Oct 2026. Truth file: 05_Exam_sealed/truth_sealed.json (SHA-256 in truth_sealed.sha256). Script: 06_Harness_scoring/score_finallaw_exam.py (reproduces every number below).

## Result in one paragraph

The law covers nine of the twelve cases. Of these, four are in the unaligned monodisperse, thick-stack or log-normal families (F1 to F4) and agree with the simulated truth to within 0.4 per cent. The two ripple-ratio cases (F10 and F11) are within 1.9 and 0.8 per cent. The three aligned-domain cases (F8, F9, F12) agree to 1.7, 1.6 and 4.2 per cent only after the order parameter and domain size had been mapped to the law's operational definitions, which was done after the truth had been seen. The mean absolute error over the nine covered cases is 1.3 per cent. The three cases outside the law (bimodal radii F7, and thickness-radius coupling F5 and F6) are not predicted by the law at all.

## Case by case

For F1 (aspect ratio 50, monodisperse) the law gives an occupied fraction of 0.02672 against a truth of 0.02663 (plus 0.3 per cent). F2 (aspect ratio 8): 0.11999 against 0.11950 (plus 0.4). F3 (aspect ratio 60, four layers): 0.07593 against 0.07625 (minus 0.4). F4 (aspect ratio 80, log-normal sigma 0.4): 0.01336 against 0.01337 (0.0). F8 (aspect ratio 40, S_dir 0.6): 0.03812 against 0.03747 (plus 1.7). F9 (aspect ratio 60, three layers, S_dir 0.6): 0.06881 against 0.06773 (plus 1.6). F12 (aspect ratio 60, three layers, sigma 0.3, S_dir 0.5): 0.05864 against 0.05627 (plus 4.2). F10 and F11 (ripple, slope 1.2, wavelength 1.0): the law's rippled-to-flat threshold ratio at fixed sheet (material) area is 1.2106 against truths of 1.188 (no polydispersity) and 1.2013 (sigma 0.4); the law has no dependence on sigma here, which the truth values suggest is a good approximation within 1 per cent.

## What this result is and is not

It is a score of the final law, with its constants fixed before the exam was simulated, on cases whose geometry the law was built for. The truth values carry roughly 0.5 to 1 per cent sampling error (24 realisations at L = 32) plus the uncertainty in the 1.015 factor, so agreement at the 0.4 per cent level on F1 to F4 is at or inside the noise. It is not independent evidence of physical correctness, because the law was calibrated to the same simulator.

It is not a blind prediction for the aligned cases. With the naive mapping (the prompt's order parameter used directly and domain size 3) the law over-predicts F8, F9 and F12 by 44.5, 37.8 and 35.3 per cent. The correct operational mapping (pair order within two radii, measured on the generator, and a domain size of six radii) was found after seeing those errors. The three aligned results are therefore a demonstration that the law has the right functional form once its parameters are properly defined, not a forecast.

The comparison is in occupied volume fraction. The law's own output is nominal, and the conversion 1 minus exp(minus phi) was needed.

Leave-one-point-out cross-validation of each constant (verification/cv_loo.py) shows no over-fitting: held-out RMS errors of 0.8, 1.7, 3.2 and 1.8 per cent for the finite-thickness, alignment, ripple and polydispersity constants against 0.8, 1.6, 3.1 and 1.7 per cent in-sample (details in 03_Provenance_Ledger.md). This tests the stability of the constants, not the functional forms; the sealed exam is the only test of the forms, and its aligned cases were not blind.

## Comparison with the baseline arms

The baseline chats scored a mean absolute percentage error of 5.0 to 9.2 per cent over their own answered cases (Opus 5.0, 6.5, 7.3; GPT 6.6, 7.2, 9.2 for single-prompt, tournament-lite and forecast-only), with a noise floor of about 5 points, and a shared under-prediction on F4 to F6. The law's 1.3 per cent on its nine cases is not comparable on equal terms: the law was calibrated on 55 simulated points, the chats were not, and the law does not cover F5 to F7. The paper may say that a law fitted by a human-AI loop with simulation access reaches the sampling-error floor on in-scope cases, and that plain-prompt chats do not; it may not say that either vendor is better.

## Open items

A fresh sealed set generated after the freeze; the effect of refitting the alignment prefactor (0.588 against 0.625) on all nine cases; extension of the law, or a declared exclusion, for bimodal and thickness-coupled distributions.
