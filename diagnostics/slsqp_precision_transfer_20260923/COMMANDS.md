# Uniform-precision full225 operations

Prepared run: `workspaces/slsqp_precision_transfer_20260923/run_v1/READY.json`.
Each manifest pins a complete implementation snapshot. Run those snapshots;
do not overwrite a completed collection or relaunch successful chemistry.
Existing old225 and repaired-H supplement are separate immutable results.

Preparation (already performed, creates a new directory only):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
scripts/slsqp_precision_transfer.py prepare \
--inventory workspaces/slsqp_precision_expansion_20260923/TRANSFER225_INVENTORY_v1.json \
--rank-qualification workspaces/native_gfn2_rank_panel_20260923/run_v2/COMPARISON.json \
--sources workspaces/union_adaptive_20260923/transfer225_v1/shard_0/proposals/manifest.json \
workspaces/union_adaptive_20260923/transfer225_v1/shard_1/proposals/manifest.json \
workspaces/union_adaptive_20260923/transfer225_v1/shard_2/proposals/manifest.json \
workspaces/union_adaptive_20260923/transfer225_v1/shard_3/proposals/manifest.json \
--agreement diagnostics/slsqp_precision_transfer_20260923/PLAN.md \
--output workspaces/slsqp_precision_transfer_20260923/run_v1
```

The four GPU allocations use `run_gpu.sbatch` with each absolute proposal
manifest as its single argument. Each GPU job executes its finite independent
searches, collects statuses, constructs that shard's three-candidate common pool
and evaluates only missing cross-MACE. Its dependent `run_cpu.sbatch` takes the
absolute pool manifest, executes its one contained native solver task manifest,
and writes `pool/collection_final.json`, including all failed cells. Actual
submission commands, IDs, dependencies and manifest hashes are saved in
`SUBMISSIONS.json`; no independent duplicate scheduler chain should be launched.

The native solver dry-run is the existing operation:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
scripts/affordable_workflow.py dry-run \
--manifest workspaces/slsqp_precision_transfer_20260923/run_v1/shard_0/pool/solvent/shard_0/manifest.json
```

After all four final collections exist, report with the separate immutable
analysis snapshot recorded in `ANALYSIS_IMPLEMENTATION.json`:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/slsqp_precision_transfer_20260923/analysis_v1/slsqp_precision_transfer_report.py \
--selection workspaces/slsqp_precision_transfer_20260923/run_v1/SELECTION.json \
--collections workspaces/slsqp_precision_transfer_20260923/run_v1/shard_0/pool/collection_final.json \
workspaces/slsqp_precision_transfer_20260923/run_v1/shard_1/pool/collection_final.json \
workspaces/slsqp_precision_transfer_20260923/run_v1/shard_2/pool/collection_final.json \
workspaces/slsqp_precision_transfer_20260923/run_v1/shard_3/pool/collection_final.json \
--baseline workspaces/union_adaptive_20260923/transfer225_v1/COMPARISON_v1.json \
--output workspaces/slsqp_precision_transfer_20260923/COMPARISON_v1.json
```
