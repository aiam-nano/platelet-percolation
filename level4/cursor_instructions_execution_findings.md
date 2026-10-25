# Instructions execution findings
_Exported on 10/25/2026 at 20:13:40 GMT from Cursor (3.27.12)_

---

**User**

Follow INSTRUCTIONS.md in this folder exactly, from Step 0 to Step 6, without asking me questions. Do not modify code/ or docs/. Do not open any other folder on this computer. Deliver findings/FINDINGS.md.

---

**Cursor**

I'll follow the adversarial review instructions in this folder from Step 0 through Step 6, working only on copies and writing the deliverable to `findings/FINDINGS.md`.

Step 0 starts with the audit packet, the provenance ledger, and the finite-size and shape reports. I'll read those next and set up the regression copy.

The packet and the three reports are in. Next I’ll copy the code, check Python, and start the line-by-line read of the simulator.

The quick suite passed 14 of 14. Next I’ll run targeted checks of the contact tests, the cell list, winding, and the law’s constants.

The winding checks matched on configurations that do percolate, and the Kell formula does not match the arc-length ratio. I’ll record that and finish the write-up while the last script completes.

The review is in `findings/FINDINGS.md`. Commands are in `review_work/COMMANDS.txt`. `code/` and `docs/` were not modified.

Disc, cylinder and spheroid contact, the cell list, and the winding detector all matched independent checks on the cases that can be run quickly. The quick regression passed 14 of 14. The thin-plate plateaus in the packaged runs sit about 1.6% (discs), 1.7% (squares) and 2.9% (triangles) above Yi and Tawerghi. Excluded volume, neighbour search, winding against a supercell, and the sphere plateau do not produce that offset, so its cause is still open.

Five findings:

- **F1.** The offset is in the data, and the law’s thin-disc limit of 0.9614 is built in: `c0` is that number and every fitted threshold was divided by 1.015. That circularity affects absolute claims only.
- **F2.** `Kell(a)` in `finallaw3.py` is not the arc-length ratio. At slope 1.2 it is 1.458 against an arc-length ratio of 1.296. In-sample ripple fits used this expression, so exam cases F10 and F11 are unchanged.
- **F3.** Polygon contact misses coplanar overlaps, and misses a vertex or edge in the other plane when the polygon lies on the negative side of the normal. Random-orientation plateaus are unaffected.
- **F4.** A pair that overlaps through two periodic images is not recorded as a wrap. Reported boxes, including exam cases F5 and F6, stay clear of that regime by about 4 length units.
- **F5.** The sphere-control sentence that quotes an extrapolation of 0.3386 against 0.3418 overstates the disagreement. The plateau of those runs is 0.3422 ± 0.0013.

Exam scores match the documents (mean absolute error 1.271% for the first law and 1.085% for the revised one). The gap tests are 22 of 24, with B2 at +18.7% and H1 at +23.6% as stated. The polydispersity matrix matches the Isihara kernel at full thickness.

---

**User**

<timestamp>Sunday, Oct 25, 2026, 8:09 PM (UTC)</timestamp>
<system_notification>
The following task has finished. If you were already aware, ignore this notification and do not restate prior responses.

<task>
kind: shell
status: success
task_id: 461893
title: Run quick regression suite
tool_call_id: call-1da9685f-0e42-473c-8316-e043b98f3cbf-18
fc_dc78cf90-c102-92bd-b087-2e448c9431eb_0
output_path: C:\Users\mbassav2\.cursor\projects\c-Users-mbassav2-Work-Work-Current-Research-UoM-AI-Google-co-scientist-Percolation-Pilot-04-Simulation-and-analysis-review-cursor\terminals\461893.txt
</task>
</system_notification>

<timestamp>Sunday, Oct 25, 2026, 8:09 PM (UTC)</timestamp>
<system_notification>
The following task has finished. If you were already aware, ignore this notification and do not restate prior responses.

<task>
kind: shell
status: success
task_id: 461894
title: Run law, exam, and document checks
tool_call_id: call-62a8c8a6-84b2-4854-b4be-31fb117d2837-69
fc_412b61a2-8f60-9984-a4e2-654c378881fd_1
output_path: C:\Users\mbassav2\.cursor\projects\c-Users-mbassav2-Work-Work-Current-Research-UoM-AI-Google-co-scientist-Percolation-Pilot-04-Simulation-and-analysis-review-cursor\terminals\461894.txt
</task>
</system_notification>

<timestamp>Sunday, Oct 25, 2026, 8:09 PM (UTC)</timestamp>
<system_notification>
The following task has finished. If you were already aware, ignore this notification and do not restate prior responses.

<task>
kind: shell
status: success
task_id: 461895
title: Run geometry adversarial tests
tool_call_id: call-62a8c8a6-84b2-4854-b4be-31fb117d2837-70
fc_412b61a2-8f60-9984-a4e2-654c378881fd_2
output_path: C:\Users\mbassav2\.cursor\projects\c-Users-mbassav2-Work-Work-Current-Research-UoM-AI-Google-co-scientist-Percolation-Pilot-04-Simulation-and-analysis-review-cursor\terminals\461895.txt
</task>
</system_notification>

<timestamp>Sunday, Oct 25, 2026, 8:09 PM (UTC)</timestamp>
<system_notification>
The following task has finished. If you were already aware, ignore this notification and do not restate prior responses.

<task>
kind: shell
status: error
task_id: 461896
title: Run wrap, face, and FSS checks
tool_call_id: call-60882d77-27a3-4f1a-8dfe-42387dbe424e-82
fc_633d779f-f957-9be2-a7a6-205743f22929_0
detail: exit_code=1
output_path: C:\Users\mbassav2\.cursor\projects\c-Users-mbassav2-Work-Work-Current-Research-UoM-AI-Google-co-scientist-Percolation-Pilot-04-Simulation-and-analysis-review-cursor\terminals\461896.txt
</task>
</system_notification>
<user_query>Briefly inform the user about the task result and perform any follow-up actions (if needed). If there's no follow-ups needed, don't explicitly say that.</user_query>

---

**Cursor**

The delayed job notifications match the review already written to `findings/FINDINGS.md`. The quick suite, the law and exam checks, and the geometry tests all finished successfully. The one failed run was the first follow-up script, which stopped because the large-box files are bare lists; the corrected script then completed and those numbers are in findings F1, F4 and F5.

