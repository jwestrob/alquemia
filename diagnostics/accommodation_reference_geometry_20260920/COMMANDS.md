# Geometry-only replay

From the repository root, with a fresh output destination:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/accommodation_reference_geometry.py scan --preparation workspaces/accommodation_goal_20260920/folds_v1/preparation_reconciled_v1.json --warning-reference /groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/xoxf_all/geometry_quality/REPORT.md --output workspaces/accommodation_reference_geometry_20260920/inventory_replay_v2.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/accommodation_reference_geometry.py mappings --inventory workspaces/accommodation_reference_geometry_20260920/inventory_replay_v2.json --output workspaces/accommodation_reference_geometry_20260920/mapping_replay_v2 --workers 1
```

These operations inspect immutable source coordinates and verify algebraic
physical maps. They call no molecular energy model. All 62 checks are already
complete under `mappings_v1`; replay is optional. The recorded batch command uses
an allocated64-CPU job to perform the same checks in parallel. No threshold,
role rule or source-metal selection is inferred from classification scores.
