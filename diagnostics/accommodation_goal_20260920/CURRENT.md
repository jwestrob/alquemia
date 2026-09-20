# Overnight accuracy research checkpoint — 2026-09-20

## Why these experiments

The release now makes source-backed PQQ scoring practical. The active goal goes
further: stop a single imperfect model from controlling a metal-specificity call,
and identify a physical adjustment that improves model reliability. The saved
reference folds test structural robustness against known classes. Whole-donor
rotations test whether diagnosed compression can be relieved with a physically
credible metal-dependent response. Unknown PLM predictions remain unlabeled.

## Current execution

- Fast-PQQ release complete: commit8dd81ba, report in pqq_fast_release_20260920.
- Matched source inventory complete: commit8ace781. All250 existingAF3 structures
  pinned, including125Ca-conditioned and125La-conditioned, no new folds.
- Root source preparation1203744 complete:233supported,17unsupported after retaining
  seven machine-last-bit replay differences as numerically identical. Only the
  comparisons are reconciled; no coordinates changed. Original failures persist.
- All834fresh native MACE endpoints complete in76.214s in a warm calculator;
  98canonical endpoints reused. The2fresh execution verification endpoints exactly
  match their archive. No scientific MACE failure; earlier1203747 failed importing
  CPU-only gemmi before any scientific evaluation and was repaired with lazy imports.
- Solvent tasks:1668new +196exact-compatible reused nativeGFN2 endpoints.
  Four disjoint case manifests run as1203797/1203798/1203799/1203800 on existing64CPU runners.
  Initial split failed path-containment validation before execution; version2
  stages actual input copies within each manifest directory. No solver results
  were generated or overwritten in either preparation failure.
- Separate physical profiles:64MACE endpoints and126/128successfulGFN2 outputs.
  Relieving compressed extraAsp donors preferentially lowersLa energy in both
  diagnosed PLM examples; this is not biological accuracy evidence. Most of the
  response is in OMOL, with the solvent term moderating it. NativeDFT1203771 now
  checks the16predeclared endpoints. Keep the failed1H4IGFN2 point unavailable.
- CC1202429 remains independent; preserve its monitor and outputs.

## Numerical replay and pooling rules

Source/PQQ/charge/water selection stays fixed. The seven canonical context
comparisons differed only by2.22e-16Å in constructed cap coordinates on another
CPU, with identical normalizedPDB bytes. The reconciliation accepts <=1e-12Å,
far below source3-decimal/core6-decimal precision, before inspecting scores.
It records raw comparisons and original statuses and performs no coordinate edit.
All other preparation failures remain unsupported.

Original25canonical geometries are replay controls. Primary pools are the other
fourLa-conditioned and fiveCa-conditioned samples per protein. The declared
medians/equal-weight mean require all members; missing preparations do not become
favorable partial ensembles. Core/composition invariants and changing context
membership remain separate. Frozen old bands are a transfer test, never refitted.

## Current files

Products: workspaces/accommodation_goal_20260920/folds_v1/.
Sourcepreparation/preparation.json is raw; preparation_reconciled_v1.json retains
last-bit acceptance; mace_v3/manifest.json is actual successful execution;
INVENTORY.json is its233-case native endpoint inventory. mace/ and mace_v2/
contain unexecuted or import-failed preparation attempts, not successful science.
The unexecuted full solvent manifest is partitioned only for throughput; use
solvent_shards_v2/ after validation and submission, not the full parent manifest.

Three actual-source adapter/reconciliation tests pass. Additional report tests
are owned by second_shell. No production rescore, PLM output change, new threshold,
mechanical correction, or ensemble promotion has occurred.

## Contained continuations

Geometry-only mapping41851ec found62flagged reference folds, but only one has
the same extraAsp compression as the two PLM examples. The others are anchorGlu
warnings and do not inherit the extraAsp validation. The sole matched-phenotype
reference (MMOL1770,Ca-conditioned sample1) is tested at fixed±0.2rad with4newMACE
and8newGFN2 calls, q0reused, in accommodation_reference_profile_20260920. No newDFT.

A separate baseline comparison is in preparation by khoury_benchmark: exact old
DFTcore recipe on all208supported noncanonical fold samples (416endpoints), with
25canonical outputs reused and17unsupported retained. This is a direct comparator
for accuracy/conditioning robustness, parallel to candidate development, with no
threshold fitting or baseline change. Confirm actualsubmission/status in that
new experiment record before running anything.

The H-only minimizer diagnosticb62b011 converged under the same generic potential
after monotonic preconditioning. It fixes silent numerical failure; that potential
still produces longHbonds, also present in its nominally successful controls.
No Hmodel change is promoted. See accommodation_controls HYDROGEN_DIAGNOSTIC_REPORT.
