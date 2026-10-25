# Review notes

Ambiguities and sources. Nothing outside this folder was opened.

## Readings I adopted

- Aspect ratio A is diameter over full thickness. That is the reading that makes phi/eta = 3/(2A), as the ledger states, and it makes Ae = A/N the aspect ratio of one cylinder of thickness N times 2/A. The moment matrix in mylaw.py passes t = 1/A into coefficients that then match the Isihara kernel for full thickness T = 2/A (see test_law_docs.py). I take that as intentional, not as a half-thickness bug.
- The thin-limit calibration targets in fitlaw.py (the Br column, about 1.01 to 1.14) are the residual multiplier (1 + kappa sigma^2), not the full phi ratio. The loss fits that residual, and the law multiplies it by the moment-matrix factor. The product matches the finite-aspect calibration points and the G1 simulations. I did not find a separate log of those five thin runs, so this reading is inferred from the loss and from G1.
- "Four times the largest bounding radius" is a sufficient condition. The condition that actually matters is that no pair can touch through two periodic images, i.e. L/2 greater than the sum of their bounding radii.
- The quick suite was run with python3, which on this machine is Python 3.14.3. The same numpy, scipy and numba are installed. A second interpreter, python, is 3.14.4. The reproducing scripts were run with python.

## Sources

- Y. B. Yi and E. Tawerghi, Phys. Rev. E 79, 041134 (2009). I used the abstract (APS and PubMed, 25 Oct 2026): thresholds 0.9614, 0.8647 and 0.7295 for circles, squares and equilateral triangles, obtained by extrapolating to zero radius. The methods section is paywalled. I did not obtain it, so I cannot certify span_open.py against their code or their equation 19, only against the author's description and that abstract.
- E. J. Garboczi, K. A. Snyder and J. F. Douglas, Phys. Rev. E 52, 819 (1995), as quoted in the Wikipedia "Percolation threshold" table (retrieved 25 Oct 2026). Oblate ellipsoids, aspect ratio 1000: occupied fraction 0.001275, reduced density nV = 0.001276. Multiplying by the aspect ratio converts that to (4/3) pi n R^3 = 1.276, which is the author's 1.27.
- Sphere value 0.3418: Lorenz and Ziff, as cited by the author; the same Wikipedia table lists 0.341889(3) and 0.289573(2) for the occupied fraction. I did not retune anything to these numbers.
- Isihara mean excluded volume of a convex body, V_ex = 2V + S M / (2 pi), used as the cylinder and spheroid check. For flat convex plates, V_ex = A P / 2.
- scipy.special.ellipe(m) is the complete elliptic integral of the second kind, integral of sqrt(1 - m sin^2) from 0 to pi/2. The arc-length ratio of z = a sin(theta) is (2/pi) sqrt(1+a^2) E(a^2/(1+a^2)).

## What I did not do

- No file outside this folder was opened or listed. The independent reimplementation was not sought.
- Long percolation campaigns were not rerun. Plateaus were recomputed from the packaged json and jsonl.
- The experimental paper texts are not in this folder. Only the law arithmetic against experimental_comparison.csv was checked.
