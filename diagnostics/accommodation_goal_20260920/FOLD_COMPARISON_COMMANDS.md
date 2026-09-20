# Frozen fold comparison operations

From the repository root, these operations parse existing outputs only. No solver,
model, geometry generation or classification fit is launched.

```bash
export FOLD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export FOLD_RUN="$PWD/workspaces/accommodation_goal_20260920/folds_v1"
```

The completed comparison is `comparison_v1.json`; its compact report is in
`diagnostics/accommodation_goal_20260920/FOLD_COMPARISON_REPORT.md`. To reproduce
collection into a fresh file:

```bash
"$FOLD_PY" scripts/accommodation_folds_compare.py compare \
 --sources diagnostics/accommodation_controls_20260920/PQQ_ALL250_SOURCES.json \
 --preparation "$FOLD_RUN/preparation_reconciled_v1.json" \
 --inventory "$FOLD_RUN/INVENTORY.json" \
 --collection "$FOLD_RUN/solvent_shards_v2/shard_0/collection_after_failure_1203797.json" \
 --collection "$FOLD_RUN/solvent_shards_v2/shard_1/collection_1203798.json" \
 --collection "$FOLD_RUN/solvent_shards_v2/shard_2/collection_1203799.json" \
 --collection "$FOLD_RUN/solvent_shards_v2/shard_3/collection_after_failure_1203800.json" \
 --reference workspaces/compact_solvation_20260920/full_v1/comparison_v1.json \
 --output "$FOLD_RUN/comparison_replay.json"
"$FOLD_PY" scripts/accommodation_folds_compare.py report \
 --result "$FOLD_RUN/comparison_replay.json" \
 --output "$FOLD_RUN/COMPARISON_REPLAY.md"
```

The two `collection_after_failure` files collect preserved completed attempts
following job-level failure. They do not rerun chemistry and preserve the two
nonconverged endpoints. Incomplete pools stay unavailable. Existing files are
immutable; use a new output name for another parser replay.
