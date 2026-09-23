# Four-call exact geometry/rank diagnosis

Completed job1211105; do not duplicate molecular work. Actual commands and
input pins are in `SUBMISSION.json` and `run_v1/READY.json`. New directories and
collection paths are write-once.

Read-only validation:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/precision_rank_geometry_20260923/run_v1/implementation/precision_rank_geometry.py \
validate --manifest workspaces/precision_rank_geometry_20260923/run_v1/rank_1/manifest.json
```

Optional read-only collection replay to a new output:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/precision_rank_geometry_20260923/run_v1/implementation/precision_rank_geometry.py \
collect --ready workspaces/precision_rank_geometry_20260923/run_v1/READY.json \
--output workspaces/precision_rank_geometry_20260923/run_v1/COLLECTION_replay_v1.json
```

The wrapper runs the existing executor concurrently on separate rank1 and rank8
manifests, then collects once. It does not rerun failed inputs, introduce seeds,
change geometry, modify thresholds, or reconstruct a classifier pool.
