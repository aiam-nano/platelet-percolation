# INSTRUCTIONS: independent reimplementation of a percolation simulator

You are the independent implementer in a verification exercise. A different AI system has already written a simulator for the model described in `SPEC.md` and wants to know whether it is correct. Your job is to write your own simulator from the specification alone, run it on the cases in `cases.json` and the test inputs in `inputs/`, and report what you obtain, whatever it is. A result that differs from what you expect, from published numbers, or from what you think the other system obtained is valuable. A result that has been nudged to look right is worthless.

Read this whole file first, then `SPEC.md`, then `cases.json`. Work autonomously from start to finish without asking questions; if something is ambiguous, pick the most natural reading, record it in `work/NOTES.md` and continue.

## Rules of independence (hard rules)

1. Stay inside this folder. Do not open, list, search or read any file outside it, including parent folders, sibling folders, the user's home directory or other projects. If your tools show you other folders, ignore them.
2. Do not search the web or any repository for code that implements this particular model, for source code named `sim2`, `first_perc`, `disc_contact`, `gjk_cyl`, or for the project that this exercise belongs to. General reference material is allowed and encouraged: textbooks and papers on geometry of intersecting discs and cylinders, the GJK algorithm, separating-axis tests, union-find with periodic offsets. Record every URL or source you consult in `work/NOTES.md`.
3. Never tune, rescale, filter or re-run anything to make a number match a published value or an expected value. If you know a published threshold for thin discs or spheres, do not use it as a target. If your first implementation disagrees with an analytic check, find and fix the bug, and log the bug and the fix in `work/NOTES.md`; that is different from tuning a result.
4. Do not edit `SPEC.md`, `cases.json`, anything in `inputs/`, or this file. All your files go in `work/`, plus the single output file `results.json` in this folder's root.
5. Do not invent data. Every number in `results.json` must come from running your code. If a case cannot be run in reasonable time, run it with fewer realisations, say so in the report and in the `nreal` field, and move on.

## Steps

Step 0. Set up. Python 3 with numpy is assumed. For speed use numba, or write the core in C, C++ or Rust and call it from Python or a script. A pure-Python inner loop will be too slow for the larger cases. Put everything in `work/`.

Step 1. Implement the two contact functions of SPEC section 2 and test them yourself. Allowed checks that need no reference values: (a) compare your thin-disc test with a brute-force test (sample many points on the circumference of each disc and test whether any lies in the other disc, and use a dense grid on one disc and the plane equation of the other); (b) compare your cylinder test with your thin-disc test as the thickness tends to zero; (c) estimate by Monte Carlo the mean excluded volume of two unit thin discs with isotropic orientations (it is known analytically to be pi squared) and of two identical cylinders with isotropic orientations (it equals 2V + S M / (2 pi), with V the volume, S the surface area, and M = pi t + pi squared R the integral of the mean curvature of the cylinder); (d) check symmetry (the result must not depend on which particle is called first). For the square and triangular plates and the spheroids of SPEC section 1: (e) estimate by Monte Carlo the mean excluded volume of two identical plates with isotropic orientations, which for flat convex plates equals A P / 2 with A the area and P the perimeter; (f) check that a spheroid with c = R gives the contact test of two spheres, and that your spheroid test approaches your thin-disc test as c tends to zero; (g) check symmetry for these types too.

Step 2. Implement the insertion process, the neighbour search and the wrapping detection of SPEC section 3, including the four recorded values N_any, N_x, N_y, N_z. Test the neighbour search against an all-pairs search on random configurations in boxes with more than three cells per side, and test the wrapping detector against a brute-force method on small boxes (for example by building the contact graph in a 3 by 3 by 3 supercell and checking whether a cluster contains two images of the same particle). Record these tests and their outcomes in `work/NOTES.md`.

Step 3. Produce the deterministic outputs. Read `inputs/pairs_inputs.json` and write, for each of the 4000 pairs, whether the thin discs intersect and whether the cylinders intersect. Read `inputs/configs_inputs.json` and write, for each of the three sets of 20 configurations, the N_any value (or -1). Then read `inputs/shapes_pairs_inputs.json` and write, for each of the three shapes (`square`, `triangle`, `spheroid`) and each of the 4000 pairs, whether the two particles intersect, and read `inputs/shapes_configs_inputs.json` and write the N_any value (or -1) for each of the 20 configurations of `sq_L8`, `tri_L8` and `sph_L8`.

Step 4. Run every case in `cases.json`. Use independent random numbers for each realisation and record the seed or generator you used. For the thin-disc cases keep the raw per-realisation arrays N_any, N_x, N_y, N_z. Compute the mean, standard error and number of realisations as in SPEC section 4. Report any realisation that failed to wrap and how many there were.

Step 5. Write `results.json` in this folder's root in exactly this structure (a template is in `results_template.json`):

```
{
  "implementation": {"language": "...", "rng": "...", "seeds": "...", "runtime_seconds": 0, "notes": "..."},
  "pairs":   {"disc": [4000 booleans], "cyl": [4000 booleans], "square": [4000 booleans], "triangle": [4000 booleans], "spheroid": [4000 booleans]},
  "configs": {"m0_t0.0": [20 ints], "m1_t0.25": [20 ints], "m1_t0.05": [20 ints], "sq_L8": [20 ints], "tri_L8": [20 ints], "sph_L8": [20 ints]},
  "stats":   {"<case name>": {"mean": 0.0, "se": 0.0, "nreal": 0, "discarded": 0}, ...},
  "raw":     {"disc": {"N_any": [...], "N_x": [...], "N_y": [...], "N_z": [...]}, "disc_L32": {...}, "disc_L48": {...}, "disc_L64": {...}, "square_L32": {...}, "square_L48": {...}, "triangle_L32": {...}, "triangle_L48": {...}}
}
```

Then run `python3 validate_results.py results.json`, which only checks the format and counts and does not compare with anything. Fix format errors; do not change numbers.

Step 6. Write `work/REPORT.md` in plain prose (no bullet points) covering: what you implemented and how; every ambiguity and the choice you made; the self-tests you ran and their outcomes, including bugs you found in your own code; the cases you ran and any deviation from `cases.json`; the numbers in a short table of case, mean, standard error, number of realisations; anything surprising; and an honest statement of how confident you are that your implementation is correct and why. Also save the full list of commands you ran in `work/COMMANDS.txt`.

Step 7. Stop. Do not look for the reference results, do not try to guess them, and do not compare your numbers to anything. The operator will run the comparison.

## What the operator will do afterwards

Copy `results.json`, `work/` and your chat transcript (exported from Cursor) back to the person who asked for this exercise. They will compare your results with the other system's results using a harness you do not have, and the code you wrote will be examined if the two disagree.

## Expected effort

All cases together are modest for compiled code: tens of minutes on a laptop for the thin-disc cases at L up to 64, a similar time for the plate cases, and a few minutes per cylinder or spheroid case. If your implementation needs hours, something is probably inefficient (check the neighbour search); do not reduce the number of realisations without recording it.
