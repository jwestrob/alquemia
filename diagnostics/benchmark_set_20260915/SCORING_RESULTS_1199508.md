# Prepared benchmark additions: completed, mixed scientific result

**All ten endpoints completed normally, with zero failures or retries.** Hans-LanM ranks strongly above the previously scored Ca-favoring GGR. Both alpha-lactalbumin preparations rank below GGR despite condition-qualified literature evidence favoring La. This batch does not establish broad La/Ca affinity discrimination.

Job **1199508**: **321 s (5m21s), 64 allocated CPUs, 20,544 allocated core-seconds**, reported CPU time 16,579 s, peak batch RSS 7,275,896 KiB, no GPU. The workflow timer separately measured 320.237676 s / 20,495.211281 allocated core-seconds. Scheduler accounting includes job startup/finish overhead. [Accounting](ACCOUNTING_1199508.txt).

## Actual results

S is in kcal/mol on the existing symmetric CN8 reporting gauge. Larger values are more La-like. **This generic v3 protocol has no calibrated absolute decision bands; negative S is not an established Ca classification.**

| Protein/site or source geometry | S | Difference from prior repaired GGR |
|---|---:|---:|
| Hans-LanM EF1 | 46.298417 | +43.240944 |
| Hans-LanM EF2 | 45.251855 | +42.194383 |
| Hans-LanM EF3 | 50.923545 | +47.866072 |
| Alpha-lactalbumin, 1F6S Ca-conditioned | −11.501773 | −14.559245 |
| Alpha-lactalbumin, 6IP9 La-conditioned | −14.485232 | −17.542705 |
| GGR, prior repaired result, reused | 3.057473 | 0 |

Hans is one ordered EF1→EF2→EF3 vector, not three independent affinity labels. Alpha's two structures are one biological observation with different source geometries and water inventories (two versus three waters). Both new targets have La-directed supporting evidence; this batch is not a balanced independent classification set.

The comparison to GGR uses the **same** `generic_peptide_amide_vertical_native_r2scan3c_v3` protocol and reporting reference. Differences are calculated directly as R_site − R_GGR, so the common aquo offset cancels. Alpha's ordering conflict therefore cannot be resolved by shifting the common reference or treating zero as a new threshold. No threshold was fitted.

Hans provides useful support for this particular protein-level comparison. Alpha is a meaningful challenge to generalization, qualified by incomplete original assay conditions and assay-to-structure construct correspondence. Its two source structures both show the conflict. These results do not isolate a cause such as solvation, geometry, protonation, water count, omitted protein response or experimental-state mismatch. No preparation, label, site selection or method was changed after inspection.

## Method and validation

Executed the exact ten prepared inputs approved in [EXECUTION_AGREEMENT.md](EXECUTION_AGREEMENT.md): baseline native r2SCAN-3c/CPCM(Water), ORCA **6.1.1**, DefGrid3, existing peptide-amide v3 preparation. Four concurrent endpoints, 16 MPI ranks each, on node-64-768g-5. This uses the baseline electronic method with the documented chemistry repair; it does not silently relabel original generic v2 results. Baseline defaults and all frozen experiments remain unchanged. No new scientific protocol ID.

All ten execution receipts, source input/XYZ/manifest hashes, executable/version, normal termination and SCF convergence were checked. Recollection exactly reproduced the saved result. Direct extraction from every archived output independently reproduced the reported endpoint energies and score algebra:

R = E_Ca − E_La; S = (R − R_aquo) × 627.509474, with energies in Hartree and R_aquo = −646.0775458314704 Hartree. Conversion occurs once. Unrounded energies and audit values are retained in [SCORING_RESULT_1199508.json](SCORING_RESULT_1199508.json).

No additional folds, environmental corrections, optimizations, threshold calibration or follow-on scientific runs were launched. The initial construction release and its unrun snapshot remain intact. The new scored ledger attaches all five scores to the two grouped evidence records and updates score exposure; all other candidates and unresolved mappings retain their previous status.

## Artifacts and next step

- Execution/collection and outputs: `workspaces/benchmark_set_20260915/ready_tasks_v4/`.
- Updated ledger: `workspaces/benchmark_set_20260915/scored_release_1199508/benchmark_manifest.json` and `.tsv`.
- Exact receipt/artifact pins: [SCORING_RESULT_1199508.json](SCORING_RESULT_1199508.json).
- Historical construction scope: [RESULTS.md](RESULTS.md).

Recommendation: **retain the baseline without claiming broad affinity validation**. Preserve alpha as a challenge case. B2/C5 and the six-blade PQQ controls remain the most directly relevant next benchmark-preparation work for expanding the paper's PQQ domain; their holo preparation still needs to be specified before scoring.

Read the completed report without starting another calculation:

```bash
cat diagnostics/benchmark_set_20260915/SCORING_RESULTS_1199508.md
```
