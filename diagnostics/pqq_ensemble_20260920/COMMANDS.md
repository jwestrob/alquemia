# Opt-in three-fold PQQ operation

Run from the repository root. The normal `standard` scorer remains the default;
`standard ensemble` is an explicit aggregation-only operation. It executes no
molecular calls and never widens the existing bands.

## Re-read the successful real run

This writes a new summary from existing completed standard outputs:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 scripts/affordable_workflow.py standard ensemble \
 --plan workspaces/pqq_ensemble_20260920/pilot_v1/plan.json \
 --collection workspaces/pqq_ensemble_20260920/pilot_v1/result_1204055.json \
 --output workspaces/pqq_ensemble_20260920/ensemble_replay_v1.json
```

Original successful output:
`workspaces/pqq_ensemble_20260920/pilot_v1/ensemble_1204055.json`.
Every output path must be new; existing immutable summaries are never overwritten.

## Declare membership before a future execution

Use the existing `fixed_core_PQQ_source_to_complete_context_v1` source schema and
released configuration. The complete executable example is [REQUEST.json](REQUEST.json).
In addition to normal case selectors/configuration, each case needs:

- `root_case_id`, `biological_group`, `source_conditioning_metal: "La"`;
- `raw_source_metal` with exact source `chain`, `resnum`, `icode`, `resname`,
  `atom` and `element: "La"`;
- its original `source_structure` path/SHA. Three distinct source hashes are
  required. Actual sequence/numbering, site roles and assembly must match.

Top-level `ensemble` contains:

```json
{
  "protocol_id": "PQQ_declared_three_La_fold_median_development_v1",
  "groups": [{
    "protein_id": "q9z4j7-pqq-la_model",
    "biological_group": "q9z4j7",
    "source_conditioning_metal": "La",
    "selection_rule": "lexicographically first three noncanonical La source IDs",
    "members": [
      "q9z4j7-pqq-la_model__conditioned_La__seed-1_sample-1",
      "q9z4j7-pqq-la_model__conditioned_La__seed-1_sample-2",
      "q9z4j7-pqq-la_model__conditioned_La__seed-1_sample-3"
    ]
  }]
}
```

All request cases must belong to exactly one group. Labels are unnecessary for
new predictions; retain `label_scope: "PQQ_prediction"`. Membership lives in the
request that `standard prepare` pins, so it cannot be changed after execution.

The following is a fresh rerun example, not a recommendation to repeat the
completed test. It requests six MACE and twelve GFN2 endpoints under existing
allocation conventions; do not submit unless that new execution is intended.
No source-code editing is required.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 scripts/affordable_workflow.py standard prepare \
 --request diagnostics/pqq_ensemble_20260920/REQUEST.json \
 --output workspaces/pqq_ensemble_20260920/fresh_repeat_v1 --mode fast-pqq

# Validate the declared source identities before calculation. Scores are unavailable.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 scripts/affordable_workflow.py standard ensemble \
 --plan workspaces/pqq_ensemble_20260920/fresh_repeat_v1/plan.json \
 --output workspaces/pqq_ensemble_20260920/fresh_repeat_v1/pre_execution.json

sbatch --job-name=PQQ-three-fold \
 --output="$PWD/workspaces/pqq_ensemble_20260920/fresh_repeat_v1/slurm_%j.out" \
 --error="$PWD/workspaces/pqq_ensemble_20260920/fresh_repeat_v1/slurm_%j.err" \
 diagnostics/pqq_ensemble_20260920/run.sbatch \
 "$PWD/workspaces/pqq_ensemble_20260920/fresh_repeat_v1/plan.json"
```

The wrapper uses the existing H200 layout:32 tasks, one CPU each, one GPU and
200000MiB host RAM. It calls unchanged standard execute/collect, then ensemble.
The source request and implementations are frozen in the plan. Explicit
`standard --mode dft-reference` is unchanged but its results are not accepted by
this context-composite ensemble operation.

If source preparation prevents a standard collection, call `standard ensemble`
without `--collection`. It retains preparation errors and returns an unavailable
median. Do not remove the failed source or select a replacement after seeing its
score. A partial component collection also yields unavailable aggregate values.
The existing standard preparer fails a multi-case scoring request when a source
is unsupported; this operation does not change that execution policy.
