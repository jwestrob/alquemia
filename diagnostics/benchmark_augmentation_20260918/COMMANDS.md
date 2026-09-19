# Operations

All operations below are local preparation/read-only verification; no energy job
is launched. Output directories/files must be new. Existing results remain pinned.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
BENCH_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
BENCH_WORK="$PWD/workspaces/benchmark_augmentation_20260918"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

"$BENCH_PY" diagnostics/benchmark_augmentation_20260918/verify_release.py \
  --manifest "$BENCH_WORK/release_v1/manifest.json" \
  --output "$BENCH_WORK/verification_replay_v1.json"

# Optional exact-policy preparation replay, not needed to use the saved inputs:
"$BENCH_PY" diagnostics/benchmark_augmentation_20260918/prepare_8gy2.py selectors \
  --source "$BENCH_WORK/sources/8GY2.cif" \
  --output "$BENCH_WORK/selectors_8gy2_replay_v1"
"$BENCH_PY" diagnostics/benchmark_augmentation_20260918/prepare_8gy2.py prepare \
  --selectors "$BENCH_WORK/selectors_8gy2_replay_v1" \
  --output "$BENCH_WORK/prepared_8gy2_replay_v1"
```

Prepared runner-compatible task manifest:

`workspaces/benchmark_augmentation_20260918/prepared_8gy2_v1/01_8GY2/pmdh_fc_holdout_8gy2_carve_manifest.json`

It contains exactly two tasks and no completed endpoint outputs. The existing
`scripts/run_orca_task_manifest.py` loads it successfully; no second workflow
system, new execution defaults or automatic submission was introduced.
