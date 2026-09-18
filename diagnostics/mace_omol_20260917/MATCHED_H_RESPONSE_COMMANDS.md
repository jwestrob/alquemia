# Matched vacuum hybrid response operations

The frozen MATCHED_H_RESPONSE_PLAN.md defines the energy, actual reused centers,
finite tasks, physical checks and accuracy decisions. Production is unchanged.
Run from the repository root, always with a fresh output path.

## Initial preparation and preflight

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_hybrid_response.py prepare --preparation workspaces/mace_omol_20260917/matched_H_prepared_v1/preparation.json --static-report workspaces/mace_omol_20260917/matched_H_report_v1/result.json --ggr-gradients workspaces/mace_omol_20260917/masked_gradient_full_report_v1/result.json --agreement diagnostics/mace_omol_20260917/MATCHED_H_RESPONSE_PLAN.md --output workspaces/mace_omol_hybrid_response_20260918/prepared_replay_v1
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_hybrid_response_20260918/prepared_v1/centers/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_hybrid_response_20260918/prepared_v1/centers/manifest.json
```

Actual `prepared_v1` contains centers(12tasks) and three whole grids(72each).
Job1201385 ran centers;1201386,1201387,1201388 run ALPHA1F6S,ALPHA6IP9,GGR1GLG.
Exact sbatch argv and manifest hashes are in each directory's submission.json.
Existing run_pilot.sbatch invokes the frozen source snapshot with one A5000,
16CPU,64474MiB. Inspect actual jobs before any retry; never duplicate a live
executor. Reuse accepted completed attempts; preserve partial/failed receipts.

All twelve new center gradients completed and scalar energies reproduce their
archives exactly (maximum difference0.0modelkcal). The two normalized whole GGR
gradients were matched exactly to actual archived inputs/receipts and reused.
Core and protein centers use the same normalized-H physical preparation.

## Initial collection and assessment

Existing wrapper writes collection_job_JOB.json on exit. Read-only collection:

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_hybrid_response_20260918/prepared_v1/centers/implementation/mace_hybrid.py collect --manifest workspaces/mace_omol_hybrid_response_20260918/prepared_v1/centers/manifest.json --output workspaces/mace_omol_hybrid_response_20260918/prepared_v1/centers/recollection_v1.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_hybrid_response.py assess --prepared workspaces/mace_omol_hybrid_response_20260918/prepared_v1/preparation.json --collections workspaces/mace_omol_hybrid_response_20260918/collections_v1.json --output workspaces/mace_omol_hybrid_response_20260918/assessment_replay_v1
```

The explicit collections_v1.json indexes all four original job collections.
Assessment requires all228 actual calls, qualifies every axis's gradient/energy
agreement, and reports the complete fine/coarse matrices and fixed sphere
solutions. Predictions are not validated scores. Unsupported curvature or
refinement results cannot trigger a native validation task.

## Conditional native validation and report

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_hybrid_minimum.py prepare --assessment workspaces/mace_omol_hybrid_response_20260918/assessment_v1/result.json --output workspaces/mace_omol_hybrid_response_20260918/minimum_replay_v1
```

Only accepted fixed predictions with unchanged typed donors yield tasks: up to
8native vacuum r2SCAN-3c EnGrad and16matching full/core masked analytic MACE.
The existing run_matched_h_quantum.sbatch dispatches through the frozen
mace_omol_matched_h.py driver; the existing GPU runner is unchanged. Validation
and collection now dispatch by the new stage, leaving original matched-H
experiments intact. No old CPCM correction, new solvent or entropy is included.

The minimum module has explicit validate,validate_quantum,collect_quantum and
report commands. Report requires --prepared (minimum preparation JSON), --quantum
(quantum manifest), --mace (learned validation manifest) and a fresh --output.
It separates actual exploratory, endpoint-validated and partition-qualified
contrasts; all four alpha/GGR comparisons must pass their original0.02margin.
Actual submitted paths/receipts will be recorded after the initial assessment.
