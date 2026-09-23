# Uniform32 fixed-pool continuation

Submitted job1211168. Do not duplicate these molecular calls. The exact command,
stage1 manifest, source inventory, old reference, frozen analysis implementation
and execution agreement are pinned in SUBMISSION.json. Both stages and final
analysis are finite operations in run.sbatch. Twelve actual earlier calls are
reused; at most756 new scalar calls execute. Missing results remain unavailable.

Read-only source/manifest validation:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/precision_pool_continuation_20260923/run_v1/stage1/implementation/precision_pool_continue_run.py \
validate --manifest workspaces/precision_pool_continuation_20260923/run_v1/stage1/manifest.json
```

When both actual collections exist, replay report/calibration to new paths:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/precision_pool_continuation_20260923/run_v1/analysis_v1/precision_pool_continue_compare.py \
--stage1 workspaces/precision_pool_continuation_20260923/run_v1/stage1/collection.json \
--stage2 workspaces/precision_pool_continuation_20260923/run_v1/stage2/collection.json \
--old-reference workspaces/slsqp_precision_expansion_20260923/REFERENCE_v1.json \
--output workspaces/precision_pool_continuation_20260923/run_v1/COMPARISON_replay_v1.json \
--reference-output workspaces/precision_pool_continuation_20260923/run_v1/REFERENCE_replay_v1.json
```

Only the designated25 stable canonical sources may define the new reference;
incomplete/unstable calibration remains an unavailable reference. Old-reference
raw transfers, all three pool geometries and all32 source statuses remain
reportable separately. No new default, full225 rescore, additional restart stage
or optimization is part of these commands.
