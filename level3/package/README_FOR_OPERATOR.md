# How to use this folder in Cursor (operator notes; the agent does not need to read this)

Open this folder as the workspace in a Cursor window with nothing else open, start a new chat, select Grok 4.7 and agent mode, allow it to run terminal commands, and paste:

Follow INSTRUCTIONS.md in this folder exactly, from Step 0 to Step 7, without asking me questions. Stay inside this folder. Do not search for existing implementations of this model and do not compare with any published value. Write everything to work/ and results.json.

Expect tens of minutes to a few hours depending on the language it chooses. When it finishes, run python3 validate_results.py results.json yourself (format only), then copy results.json, work/ and the exported chat transcript out of the folder. Run the comparison from the held-back folder (python3 compare_reimplementation.py path/to/results.json) and read README_COMPARISON.md there. Do not put the held-back folder, the old specification in 08_Drafting, or the old verification/reimpl folder anywhere Cursor can see them.


Added 25 Oct 2026 (before the folder was given to Grok): SPEC.md, INSTRUCTIONS.md, cases.json, results_template.json and validate_results.py now also cover zero-thickness square plates, zero-thickness equilateral triangles and oblate spheroids of aspect ratio 1000, with new inputs in inputs/shapes_*. The previous versions are in ../level3_package_backup_24Oct. The held-back reference for the additions is in reimplementation_HELDBACK (shapes_*.npz, generate_reference_shapes.py, reference_stats.json); compare_reimplementation.py now has sections A2, B2 and C2 for them. Expect somewhat longer runs than before.
