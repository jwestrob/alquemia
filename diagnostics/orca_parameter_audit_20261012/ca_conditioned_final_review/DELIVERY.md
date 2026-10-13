# Corrected La/Ca recalibration and Ca-conditioned transfer

12 October 2026. The corrected expanded-context model retains all 25 canonical
reference decisions and transfers well to the declared Ca-conditioned structures.
Accommodation provides a modest improvement over its separately calibrated static
counterpart, with one persistent wrong call, one abstention and explicit coverage
losses. Production has not changed.

## What this establishes

The test uses the published-parameter GFN2 solvent components, the unchanged MACE
checkpoint, the existing expanded context and the same finite shared geometry
pool. Each metal chooses among the actual origin, Ca proposal and La proposal.
The operational policy retains the origin unless improvement exceeds 0.1 kcal/mol.
All decisions use the previously frozen canonical calibration for that policy;
this Ca-conditioned test did not supply new thresholds.

| Population | Static: correct / wrong / inconclusive / unavailable | Accommodated |
|---|---:|---:|
| Canonical calibration, 25 | 25 / 0 / 0 / 0 | 25 / 0 / 0 / 0 |
| Ca-conditioned sources, all 125 | 104 / 1 / 4 / 16 | 105 / 1 / 1 / 18 |
| Same 107 available sources | 102 / 1 / 4 / 0 | 105 / 1 / 1 / 0 |

Accommodation resolves A0A3F2YLY8 samples 0, 1 and 3 and MMOL_1770 sample 1.
MMOL_1770 sample 0 changes from correct to inconclusive. The known A0ACD6B9F2
sample-4 Ca-conditioned structure remains wrongly Ca-classified. These are all
seed-1 structures. The full denominator preserves 16 preparation exclusions and
two P12293 searches interrupted by the earlier node failure; their valid origins
do not silently stand in for missing accommodated scores.

The canonical operational gap is 3.472171047 kcal/mol. The older corrected
strict208 calibration still fails separation and remains a separate negative
result. Earlier La-conditioned results and 94 available correct three-source
summaries remain separate evidence: this test defines no new Ca-conditioned
triple rule. Related folds are correlated development samples, not independent
biological validation. These contrasts are not binding free energies.

## Execution and verification

MACE recovery lanes 2660/2661 completed in 3h09m18s and 3h15m31s on seven CPUs
each; collector 2662 took 13 seconds. The corrected solvent stage 2663 completed
in 2h06m44s on 12 CPUs: 1,264 xTB calls plus 1,264 independent tblite checks,
88,725.616 measured child CPU-seconds and 91,248 scheduler-allocated CPU-seconds.
These are stage costs, including the qualification backend, not single-protein
latency or total campaign cost. Reuse, earlier attempts and interruptions are
retained separately in the execution records.

The final saved-artifact review runs under Slurm and makes no molecular calls.
Attempt 2832 stopped on an excluded-source metadata key before scientific review;
its failure is preserved. Version 2 joins using existing source IDs while retaining
null prepared-case IDs. The reporting node then returned Slurm protocol errors,
retained stale running allocations and had no review worker in process inspection.
It resolved to a different hostname from its scheduler entry. The unchanged
version-2 analysis was submitted as job 2834 on registered node-64-768g-3, with
separate outputs and completion watcher 3426029. No molecular work was repeated
and no running job was cancelled. Its actual verification receipt, tables and
generated report are under
`workspaces/orca_parameter_audit_20261012/ca_conditioned_final_review_v3/`.
Job 2834 completed successfully in 54 seconds on one CPU. The actual review
reparsed all 1,300 unique inputs, checked all 650 solvent components and reproduced
every recorded selection and decision. Maximum independent component discrepancy
was 1.339572577e-7 Eh (frozen limit 1e-6), and population discrepancy was
1.049698817e-7 e (limit 1e-5). All original absolute-energy v1 discrepancies remain
failed; qualification concerns the matched solvent differences, not absolute
cross-implementation totals. The results and all 125 rows are committed alongside
this delivery in RESULT.json, REPORT.md, ALL125.tsv and EXCEPTIONS.tsv.

The retained native receipts record 3,461 MACE evaluations and 263,376.255 inference
CPU-seconds across the proposal history. Interrupted forwards without terminal
receipts retain unknown costs; the solvent-stage timing is not the full campaign.

All 58 scored Ca-class sources are correctly classified. Among 49 scored La-class
sources, 47 are correct, one is wrong and one is inconclusive. The separate earlier
La-conditioned noncalibration population contains 104 correct accommodated calls.
Do not combine overlapping populations into an independent-validation denominator.

Primary collection SHA256:
`304486558f3bf821d81dd270758bd0136b35440004c7a87baf2bc9892ce9cec4`.

## Decision

Retain the corrected expanded candidate for the next versioned release assessment.
This round establishes useful corrected reference fidelity and a limited benefit
from accommodation. It does not justify moving thresholds to hide the remaining
error, claiming uniform structural robustness, or transferring the old native
score scale. Keep the historical results and production identities intact.

The La/Ca molecular campaign is complete. Existing LanM jobs 2665/2666 and their
comparison 2669 remain monitored; no new La/Dy work was launched for this delivery.
