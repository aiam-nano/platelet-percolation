# Specification: continuum percolation of penetrable platelets in a periodic box

This document defines what must be implemented. It contains no source code and no reference results. If something is ambiguous, choose the most natural reading, write the choice down in `work/NOTES.md`, and carry on; do not try to find out what the "expected" answer is.

## 1. Model

Particles are circular cylinders (the thin-disc case is the zero-thickness limit). A particle has a centre p, a unit normal n (the cylinder axis), a radius R, and a full thickness t, so the solid is the set of points whose component along n lies in [-t/2, t/2] and whose component perpendicular to n has length at most R. Lengths are in units of the mean radius. Particles are fully penetrable: they may overlap freely. Two particles are connected if their closed solids share at least one point (direct contact, no tunnelling distance, no shell).

Three further particle types are used in some cases (added 25 Oct 2026, before the package was given to the second system). A particle of any type has a centre p and an orthonormal frame (e1, e2, n); the particle's own coordinates are measured along e1, e2 and n from p.

Square plate: a zero-thickness square in the plane through p spanned by e1 and e2, with side s equal to the square root of pi (area pi), sides parallel to e1 and e2, centred on p; its vertices are at (plus or minus s/2, plus or minus s/2) in the (e1, e2) coordinates.

Triangular plate: a zero-thickness equilateral triangle in the same plane, with side a equal to the square root of (4 pi / square root of 3) (area pi), its centroid at p, and its three vertices at distance a / square root of 3 from p, at in-plane polar angles 90, 210 and 330 degrees measured from e1 towards e2.

Oblate spheroid: the solid set of points whose component along n is x_n and whose component perpendicular to n has length x_p, with (x_p / R) squared plus (x_n / c) squared at most 1; R is the equatorial radius (R = 1 in the cases below) and c the polar semi-axis. The aspect ratio is R / c.

For these types the bounding radius is the distance from p to the farthest vertex (square plates: s / square root of 2; triangular plates: a / square root of 3) or the larger of R and c (spheroid). The requirement that L exceed four times the largest bounding radius applies as before.

The simulation box is a cube of side L with periodic boundary conditions in all three directions. Two particles interact through their nearest periodic images only (the code must require L larger than four times the largest bounding radius, where a particle's bounding radius is sqrt(R squared + (t/2) squared)).

## 2. Contact tests

Implement two functions.

Thin discs (t = 0): an exact test. Two zero-thickness discs intersect if and only if they share a point; two non-parallel planes meet in a line, and the discs intersect when the chord of that line inside the first disc overlaps the chord inside the second disc. Handle the parallel and nearly parallel cases explicitly (coplanar discs intersect when the distance between their centres is less than the sum of the radii; parallel non-coplanar discs do not intersect).

Cylinders (t > 0): an exact test between two convex bodies. Any exact method is acceptable (separating-axis tests with the full set of candidate axes, GJK with an analytic support function, or a distance minimisation). Do not approximate a cylinder by a polyhedron or by spheres. Tangent contacts have measure zero and are not important, but results must not depend on the tolerance in any visible way.

Square and triangular plates: an exact test. Two zero-thickness convex polygons in non-parallel planes meet in a line segment if and only if the chord of the line inside the first polygon overlaps the chord of the line inside the second; compute each chord by clipping the polygon against the plane of the other. Parallel planes: treat as not intersecting (exactly coplanar polygons occur with probability zero in every test below, and need not be handled).

Oblate spheroids: an exact test between two convex bodies. Any exact method is acceptable (GJK with the analytic support function of an ellipsoid, or an exact ellipsoid-ellipsoid separation test). Do not approximate a spheroid by a polyhedron or by spheres, and note that the cases below include a spheroid of aspect ratio 1000, so the method must stay accurate for extremely flat bodies.

## 3. Insertion process and percolation criterion

Particles are inserted one at a time, in a fixed order, at positions uniform in the box. After each insertion the new particle is tested against all earlier particles, using a neighbour search (a cell list is recommended; verify it against an all-pairs search). Connected clusters are maintained incrementally. The threshold of a realisation is the number N of particles inserted at the moment a cluster first contains a closed path of contacts whose net displacement is a non-zero vector in units of the box side, that is, the cluster wraps around the periodic box. N is the count including the particle whose insertion closed the wrapping path.

A standard way to detect wrapping is a union-find structure that stores, for every particle, the integer lattice offset of its image relative to its cluster root; when a new contact joins two particles already in the same cluster, the offset mismatch is the winding of the closed path, and a non-zero mismatch means the cluster wraps. If you use another method, say so.

Record four numbers per realisation, not just one. `N_any` is the N defined above (the first non-zero winding in any direction). `N_x`, `N_y`, `N_z` are the first insertion counts at which some closed path of contacts has a non-zero winding component along x, y or z respectively. Always `N_any = min(N_x, N_y, N_z)`.

If no wrapping occurs when all particles of a realisation have been inserted, discard that realisation and report how many were discarded (there should be none if you insert enough particles; the number to insert is given per case).

## 4. Quantities to report

For thin discs, the reduced density eta = (4 pi / 3) * N * <R cubed> / L cubed, where <R cubed> is the mean of R cubed over the N inserted particles (for equal unit discs this is simply (4 pi / 3) N / L cubed).

For square plates, triangular plates and spheroids, the reduced density eta = (4 pi / 3) * N * <r cubed> / L cubed with r = 1 for every particle (the plates have area pi, that is, the same area as a unit disc, and the spheroids have R = 1). It is therefore the same formula as for the unit thin discs: eta = (4 pi / 3) N / L cubed.

For cylinders, the nominal volume fraction phi = N * <V> / L cubed, where V = pi R squared t is the particle volume and <V> is the mean over the N inserted particles.

For each case report the mean over realisations, the standard error of the mean, and the number of realisations. Also report the raw per-realisation values of N_any, N_x, N_y, N_z for the thin-disc cases (the raw arrays are needed to compare alternative estimators).

## 5. Systems

Radius R = 1 unless stated; thickness t is the full thickness; the aspect ratio A is defined as 2R/t.

`disc`, `disc_L32`, `disc_L48`, `disc_L64`: thin discs, R = 1, t = 0, random normals uniform on the sphere.

`cyl_A8`: R = 1, t = 2/8. `cyl_A50`: R = 1, t = 2/50. Random normals uniform on the sphere.

`stack3_A60`: R = 1, t = 3 * 2/60 = 0.1 (a stack of three layers is modelled as one thicker cylinder).

`lognormal0.4_A80`: radii drawn independently from a log-normal distribution with standard deviation of ln R equal to 0.4 and mean of ln R equal to minus 0.08 (so that the mean of R is 1), samples with R above 6 discarded and redrawn by simply dropping them, t = 2/80 for every particle (independent of R). Random normals.

`aligned_k6.312_ell6_A40`: R = 1, t = 2/40. Orientations are locally correlated as follows. Place L cubed / 216 seed points uniformly at random in the box (a domain size of 6 radii). Give each seed an independent random unit vector, its director. Assign every particle to its nearest seed using the periodic minimum-image distance (a Voronoi domain). Draw each particle's normal from a Fisher distribution about the director of its domain with concentration kappa = 6.312: the cosine of the polar angle is c = 1 + ln(u + (1 - u) exp(-2 kappa)) / kappa with u uniform on [0, 1], and the azimuth is uniform. Before using the generator, verify that for kappa = 6.312 the average of P2(cos theta) about the domain director is about 0.60, where P2(x) = (3 x squared - 1)/2; put your measured value in `work/NOTES.md`.

`square_L32`, `square_L48`: zero-thickness square plates (section 1), the frame of every particle a uniformly random rotation (equivalently, n uniform on the sphere and the in-plane orientation uniform in angle). `triangle_L32`, `triangle_L48`: the same with triangular plates. `spheroid1000_L32`, `spheroid1000_L48`: oblate spheroids with R = 1 and c = 1/1000, n uniform on the sphere.

## 6. Cases to run

The cases, box sizes, numbers of realisations and numbers of particles to insert are in `cases.json`. The number to insert is an upper limit (about twice the expected threshold); stop each realisation as soon as it wraps.

## 7. Deterministic test inputs

`inputs/pairs_inputs.json` (also as `.npz`) contains 4000 pairs of particles. For each, the first particle is at the origin with normal n1, radius R1 and half-thickness h1, and the second is at displacement d with normal n2, radius R2 and half-thickness h2 (so full thickness t = 2 h). Return, for each pair, whether the two thin discs (ignoring the thicknesses) intersect, and whether the two cylinders (with the thicknesses) intersect. No periodic images are involved in these pairs.

`inputs/configs_inputs.json` (also as `.npz`) contains three sets of 20 small configurations in a box of side 6, 120 particles each, in a fixed insertion order, with positions and normals given. The sets are thin discs (t = 0), cylinders with t = 0.25 and cylinders with t = 0.05, all with R = 1. For each configuration return the N at which wrapping first occurs (N_any), or -1 if it never occurs with the 120 particles. Insert the particles in the order given.

`inputs/shapes_pairs_inputs.json` (also as `.npz`, keys `<shape>_<name>`) contains, for each of three shapes (`square`, `triangle`, `spheroid`), 4000 pairs. For pair i of a shape, particle 1 is at the origin with orthonormal frame (e1_1[i], e2_1[i], n1[i]) and particle 2 is at displacement d[i] with frame (e1_2[i], e2_2[i], n2[i]); the frames are given explicitly so that no convention about in-plane orientation is needed. The spheroid pairs also give the equatorial radii R1, R2 and the polar semi-axes c1, c2, which differ from pair to pair (aspect ratios from 1 to 1000); the plate pairs use the plates of section 1. Return, for each pair of each shape, whether the two particles intersect. No periodic images are involved.

`inputs/shapes_configs_inputs.json` (also as `.npz`, keys `<set>_s<k>_pos`, `_e1`, `_e2`, `_n`) contains three sets of 20 configurations in a box of side 8, 250 particles each in a fixed insertion order, with positions and frames given: `sq_L8` (square plates), `tri_L8` (triangular plates) and `sph_L8` (spheroids with R = 1 and c = 0.001, frame vector n only). For each configuration return the N at which wrapping first occurs (N_any), or -1 if it never occurs with the 250 particles. Insert the particles in the order given.

## 8. What is deliberately not specified

The random number generator, the data structures, the language, and the cell size of the neighbour search are yours to choose. The programme should run the whole of `cases.json` in no more than a few hours on a laptop; if a pure-Python version is too slow, use numba, C, C++ or Rust.
