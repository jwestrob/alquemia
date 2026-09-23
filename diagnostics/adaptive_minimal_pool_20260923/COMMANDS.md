# Reproduce the archive-only minimal-pool experiment

These operations make no molecular calls. Use new output paths: existing records
are deliberately immutable. Calibration must precede the transfer comparison.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
NIKASHA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python

"$NIKASHA_PY" scripts/adaptive_minimal_pool.py calibrate \
  --collection workspaces/adaptive_completion_20260922/original30_pool_v1/final_collection.json \
  --agreement diagnostics/adaptive_minimal_pool_20260923/PLAN.md \
  --output workspaces/adaptive_minimal_pool_20260923/REFERENCE_replay.json

"$NIKASHA_PY" scripts/adaptive_minimal_pool.py compare \
  --collection workspaces/adaptive_completion_20260922/primary225_pool_v1/final_collection.json \
  --reference workspaces/adaptive_minimal_pool_20260923/REFERENCE_replay.json \
  --full-comparison workspaces/adaptive_completion_20260922/primary225_pool_v1/final_comparison.json \
  --output workspaces/adaptive_minimal_pool_20260923/COMPARISON_replay.json

"$NIKASHA_PY" -m unittest discover -s tests -p test_adaptive_minimal_pool.py -v
```

The completed actual artifacts are `REFERENCE_v1.json` and `COMPARISON_v1.json`
in that workspace. The implementation is a report-only candidate-pool replay,
not a replacement production source-to-score command. The separately owned
two-source recovery adds direct preparation/execution of the smaller pool.
