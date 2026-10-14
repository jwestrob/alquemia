# Inconclusive benchmark review — 14 October 2026 UTC

Jacob requested review of the fixed-charge candidate's inconclusives. This is
saved-artifact analysis; production, bands, labels and all historical predictions
remain unchanged. No new molecular evaluations or folds are authorized by these
diagnostic completion events.

**Complete: read [the scientific report](REPORT.md).** Root score job4037 completed
all48 declared source/context comparisons. Geometry job4038 completed all48 sources
and all three saved candidates per source; corrected contact summary4041 is complete. The
nine target contexts span five proteins; four belong to the recovered Ca-source
125 and five additional contexts are La-conditioned ExaA/PedE controls.

Read `scores_v1/execution/RESULT.json` and `geometry_v1/execution/RESULT.json`
under `workspaces/pqq_inconclusive_review_20261014/`. All inputs, code hashes,
Slurm commands and actual completion-wake commands are pinned in that workspace.
The geometry result's `closest_mapped_nonbonded_heavy_pairs` field is invalid:
the kinematic graph omits fixed PQQ internal bonds. Do not use it as a clash
measurement. Other geometry/confidence/motion outputs are available. The separately
versioned `geometry_contacts_v2/execution/RESULT.json` supersedes only that field,
preserving the original output and all scientific scores.

All agent tasks and owned analysis jobs are complete. Three one-CPU jobs used53
allocated CPU-seconds, no new energies or folds. Root owns integration. The proposed
next physical test is a matched, fixed-geometry omitted-packing comparison in ExaA;
it has not been launched. No new energy calculation, optimization, refold or recalibration is queued.
Unrelated scheduler-listed reporting jobs2669/2832/2833 remain untouched.

The four current Ca-source ambiguities are ExaF two-Ca4, XoxF4-2 Ca0, ExaA
two-Ca3 and XoxF1 Ca4. Five additional ambiguous ExaA/PedE La-conditioned controls
and the five existing inconclusive triples remain in the denominator. No common
geometry/confidence/bound rule explains them. Reopening older terminal events
does not authorize repeat calculations. Next read-only command:

```bash
cat diagnostics/pqq_inconclusive_review_20261014/REPORT.md
```
