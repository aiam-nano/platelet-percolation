# Notes

No URL or external source was consulted. The contact tests, the GJK iteration and the periodic union-find were written from the specification and from standard geometry: the chord of a plane-plane line, the support function of a cylinder and of an ellipsoid, the closest point on a triangle, and a union-find that stores a lattice offset to the cluster root.

## Ambiguities and the choice made

Closed discs. Section 1 connects particles when the closed solids share a point. Section 2 says coplanar discs intersect when the centre distance is less than the sum of the radii. The implementation uses a closed test: coplanar discs intersect when the in-plane centre distance is at most R1+R2, and non-parallel discs intersect when the closed chords overlap, including a shared endpoint.

Parallel plates. Section 2 says parallel planes, including the exactly coplanar case, are treated as non-intersecting. The plate test returns no contact when sin squared of the angle between the normals is below 1e-20. A vertex is treated as lying in the other plane when its signed distance is at most 1e-12. Crossing squares that touch at one boundary point are reported as intersecting.

Nearly parallel discs. The chord test is used whenever sin squared of the angle is at least 1e-24. Below that, the planes are treated as parallel. A separation along the normal larger than 1e-10 times (1+R1+R2) is a miss; a smaller separation is treated as coplanar. For radii of order 1 this threshold sits far below any distance the Monte Carlo checks can see.

Cylinders and spheroids. Intersection is decided by GJK on the Minkowski difference, using the analytic support function. The origin is reported inside when its distance to the current simplex is at most 1e-9. A support evaluation that cannot move the closest point by more than that tolerance is reported as a miss. Bounding spheres reject pairs before GJK. The same tolerance therefore counts a gap smaller than about 1e-9 as contact. Deterministic certification, below, is how that choice was checked.

Insertion cap. The strings "0.40*L^3" and "2.5*L^3" are read as an upper bound, and the number of particles is the greatest integer not above the value. The caps used were 5529, 13107, 44236 and 104857 for L = 24, 32, 48 and 64, and 34560 for every L = 24 cylinder case.

Which N enters the reported mean. The mean is the mean of eta or of phi evaluated at N_any, with the sample standard deviation (divisor n-1) divided by sqrt(n) as the standard error. eta is (4 pi / 3) N / L^3 for thin discs, plates and spheroids. phi is N times the mean particle volume over L^3, which for the log-normal case is pi t times the sum of R squared over the N particles, divided by L^3.

Directional thresholds. For cases with record_raw the insertion continues until a wrap has been seen along x, along y and along z, or the cap is reached. N_any remains the first of those events. Every raw realisation wrapped in all three directions inside the cap, so none was extended and none was censored. Cases without raw output stop at the first wrap.

Minimum image. The image shift of component delta is s = -floor(delta/L + 0.5), so the wrapped component lies in [-L/2, L/2). A tie at +L/2 is mapped to -L/2.

Random numbers. Each realisation uses numpy Generator PCG64 seeded by SeedSequence([20261025, case_index, realisation_index]), with case_index the order in cases.json starting at 0. Thin discs, equal cylinders and spheroids draw, per particle, three uniforms for the position and then a z-cosine and an azimuth. Plates draw one further uniform spin. Log-normal particles draw the position, then normal deviates for ln R until R is at most 6, then the orientation. The aligned case draws, for each of the 64 seeds, a position and a director, then for each particle a position and the two Fisher uniforms. The Fisher cosine is the formula in the specification. The in-plane frame of a spheroid is completed with spin 0, because the solid is symmetric about its normal. Plate and Fisher frames use the cross product of the normal with the z axis when |n_z| < 0.9, and with the x axis otherwise.

Aligned seeds. The seed count is round((L/domain_size)^3) = 64, which equals L^3/216 for this case. The nearest seed uses the same minimum-image shift as the contacts.

Log-normal box. R is capped at 6 and L is 24, so L is not strictly larger than four times the bounding radius of a particle with R near 6 (the half-thickness adds a little to R). Contacts still use only the nearest image, which is what the specification requires.

Neighbour cells. The stencil radius is floor(interaction range / cell edge) + 1, with the range equal to the sum of bounding radii plus 1e-8. A cell edge a little above the typical interaction range is chosen, and at least about four cells on a side when the box allows it. If the stencil covers the whole box, each cell is visited once.

Union-find. shift[i] is the lattice image of particle i that is connected to the root's image in cell 0. A new contact between i and the image s of j produces winding s + shift[i] - shift[j] when the two particles already share a root. A nonzero winding is recorded and that edge is not unioned, so the offsets stay single-valued. The first insertion index at which a component becomes nonzero is N_x, N_y or N_z, counting from 1.

## Self-tests

Disc witnesses on 2000 random pairs all lay in both discs, the predicate was symmetric, and a 64-point circumference sample never found a point in the other disc when the predicate said miss. The Monte Carlo excluded volume of unit discs was 9.932 with standard error 0.037, against pi squared = 9.8696 (400000 samples in the cube [-2, 2]^3).

Cylinder hand cases (lateral overlap and gap, axial overlap and gap, perpendicular crossing and a separated skew pair) all matched. On 1500 random cylinders every reported hit had a witness inside both cylinders, and the predicate was symmetric. A zero-thickness cylinder agreed with the disc test on every one of 2000 pairs whose disc answer was stable under a 1e-4 centre shift. The fraction of disagreements with the disc test at half-thickness 1e-2, 1e-3 and 1e-4 was 0.005, 0 and 0 (3000 pairs each). The excluded volume at R = 1 and t = 0.5 was 20.303 with standard error 0.082, against 2V + SM/(2 pi) = 20.302 (150000 samples).

The separation-direction check in that same cylinder sample initially reported 331 failures. Those misses were bounding-sphere rejections, which return a zero vector rather than a direction. Wherever a direction was returned, the analytic gap matched the GJK distance to 2e-8. The production certification below uses that distinction and found no failures.

Plates: perpendicular squares through the origin intersect, the same squares shifted by more than a side do not, exactly coplanar overlapping squares are reported as non-intersecting, and squares that touch at one point intersect. On 800 random squares and 800 random triangles every hit had a witness in both polygons, and the predicate was symmetric. Square excluded volume was 10.981 with standard error 0.126 against AP/2 = 11.137 (80000 samples). Triangle excluded volume was 13.181 with standard error 0.194 against AP/2 = 12.693. The triangle estimate is about 2.5 standard errors high. The witness check rules out false hits, and the deterministic triangle pairs later certified cleanly, so this was left as Monte Carlo noise.

Spheroids with c = R matched the sphere test on all 2000 pairs whose centre distance was more than 1e-6 from R1+R2, and the predicate was symmetric. The fraction of disagreements with the unit-disc test at c = 0.05, 0.01 and 0.002 was 0.008, 0.0013 and 0 (1500 pairs each). The sphere-witness check in that run reported 441 failures. The boolean matched the sphere test; the witness was NaN because a zero barycentric weight was multiplied by uninitialised simplex memory. The simplex buffers were changed to zeros before the production run. After that change, all 4000 deterministic spheroid pairs certified (witness inside both ellipsoids on a hit, separating gap matching the GJK distance on a miss).

Fisher concentration 6.312. The mean of P2(cos theta) was 0.6011 over 200000 cosine samples and 0.6006 over 100000 normals drawn about a fixed director.

Wrapping. Seven coplanar unit discs at x = 0.5, 2, 3.5, 5, 6.5, 8, 9.5 in a box of side 10 wrap at insertion 7, and only in x. Twelve random disc configurations of 36 particles in a box of side 12 gave the same (N_any, N_x, N_y, N_z) on a one-cell grid, a 5-cell grid and an 8-cell grid. At the reported wrapping insertion, and at the previous insertion, those values agreed with a brute-force 3 by 3 by 3 supercell (zero disagreements).

## Deterministic certification

On the 4000 input pairs the disc, cylinder, square, triangle and spheroid predicates produced a valid witness for every hit and, for cylinders and spheroids, a separating direction consistent with the support function for every miss that was not a bounding-sphere rejection. The failure counts were all zero. The first configuration of each of the six configuration sets produced the same wrapping numbers on the production cell list and on a single-cell all-pairs list.

Intersection counts on the 4000 pairs were 495 discs, 1086 cylinders, 339 squares, 184 triangles and 835 spheroids. Every one of the 20 configurations in each of the six sets wrapped inside the supplied particle list.

## Cases

No realisation was discarded. No raw realisation needed particles beyond the cap, and no directional threshold was censored. The means are in REPORT.md and in results.json. `python validate_results.py results.json` printed FORMAT OK.
