# Explicit prepared-core response scoring

`scripts/response_probe_score.py` is the reusable opt-in scoring entry point.
It accepts a real prepared core, completed native MACE collection, compatible
native DFT task manifest, case ID, saved response evaluation/model, and baseline
release. It performs no fit or molecular calculation and never changes defaults.

The function `score(preparation, mace_collection, dft_manifest, case_id,
evaluation, baseline_release, output)` has the same explicit arguments. Output
is a new JSON file; existing files are never overwritten. The caller may use
this function after its normal preparation/DFT and new native MACE-core tasks.

```bash
OPENBLAS_NUM_THREADS=1 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/response_probe_score.py \
  --preparation workspaces/benchmark_augmentation_20260918/prepared_8gy2_v1/01_8GY2/pmdh_fc_holdout_8gy2_carve_manifest.json \
  --mace-collection workspaces/response_probe_20260919/continuation_v2/mace/collection_job_1202089.json \
  --dft-manifest workspaces/response_probe_20260919/continuation_v2/dft_8gy2_retry_v1/manifest.json \
  --case-id 8GY2 \
  --evaluation workspaces/response_probe_20260919/evaluation_v2.json \
  --baseline-release diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json \
  --output workspaces/response_probe_20260919/8GY2_scoring_replay_v1.json
```

## Output contract

- `status`: `complete` or `unavailable`, with explicit reason. CLI exits nonzero
  on unavailable; failed force calculations cannot produce a successful head.
- `baseline`: native `R_kcal_mol`, canonical reporting `S_kcal_mol`, and its own
  released decision/bands. Generic/direct cores retain raw R without a borrowed S.
- `features`: raw DFT_R plus mean radial differential donor gradient and the
  descriptive tangential fraction; all units appear in field names/measurement.
- `response_heads.DFT` and `.DFT_radial`: separate saved PQQ-model logits/classes
  at frozen zero-logit threshold. These are research functional/association
  classifications, not affinities or probabilities. Direct/generic cores export
  physical features while absolute cross-target classes remain unavailable.
- `DFT_endpoints`, `force_artifacts`, `sources`, `physical_mapping`, `donors` and
  `response_measurement`: source/receipt/coordinate/charge/mapping details.
- `new_fits=0`, `new_molecular_calls=0`, `baseline_changed=false`; relaxation and
  entropy corrections are null. A retained baseline alongside a failed head is
  not reported as head success.

Native forces must come from the same pinned float64 MACE-POLAR-medium vacuum
checkpoint/software/kernel used for the fitted features. Core positions, atom
order, charges, physical caps and explicit waters must exactly match DFT. The
canonical head accepts only the unchanged fixed-core v3 PQQ microstate/water
policy; the DFT runner's actual receipt, software version and input method are
checked. The learned head does not make an excluded heme or whole protein valid.

The completed six-task continuation uses `scripts/response_probe_run.py` as a
small manifest validator/collector around the unchanged preserved native worker
and existing execute/cache/lock implementation. The general scoring operation
does not depend on that exact task set or those three protein IDs.

For scanner integration, keep this head opt-in beside existing baseline fields.
The observed GGR structural sensitivity and any PQQ transfer failure must remain
visible; the interface being runnable is not evidence of broad validation.
