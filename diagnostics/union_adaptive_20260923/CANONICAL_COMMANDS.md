# Canonical continuation operations

Use the declared workspace; do not resubmit completed molecular jobs. These
commands validate preparation and later replay reporting into a new output.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
UNION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
UNION_RUN=$PWD/workspaces/union_adaptive_20260923
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$UNION_PY" scripts/union_adaptive.py validate \
  --manifest "$UNION_RUN/canonical_proposals_v1/manifest.json"
"$UNION_PY" scripts/nikasha_pool.py validate \
  --manifest "$UNION_RUN/canonical_pool_v1/manifest.json"
"$UNION_PY" scripts/union_adaptive.py canonical_reference \
  --collection "$UNION_RUN/canonical_pool_v1/collection_final.json" \
  --inputs "$UNION_RUN/CANONICAL_INPUTS_v1.json" \
  --agreement diagnostics/union_adaptive_20260923/CANONICAL_PLAN.md \
  --output "$UNION_RUN/CANONICAL_REFERENCE_replay.json"
```

The reference replay requires the actual terminal collection. It uses exactly25
original designated calibration sources and excludes all three crystals. Both
mathematical and operational variants retain the existing0.02kcal minimum gap.
No old band is silently inherited. Replaying only changes the freeze timestamp,
not actual input energies or membership.

The exact source selection and molecular preparation were produced with:

```bash
"$UNION_PY" scripts/union_adaptive.py canonical_inputs \
  --calibration workspaces/consistent_context_20260922/calibration28_v1/collection_final_v2.json \
  --pilot "$UNION_RUN/pool_v1/collection_final.json" \
  --agreement diagnostics/union_adaptive_20260923/CANONICAL_PLAN.md \
  --output "$UNION_RUN/CANONICAL_INPUTS_replay.json"
"$UNION_PY" scripts/union_adaptive.py prepare \
  --preparation workspaces/consistent_context_20260922/prepared_v1/PREPARATION.json \
  --crystals workspaces/consistent_context_20260922/crystal_controls_v1/CRYSTALS.json \
  --calibration workspaces/consistent_context_20260922/calibration28_v1/collection_final_v2.json \
  --transfer workspaces/consistent_context_20260922/primary225_v1/collection_final_v1.json \
  --inputs "$UNION_RUN/CANONICAL_INPUTS_replay.json" \
  --source workspaces/adaptive_completion_20260922/original30_v1/manifest.json \
  --agreement diagnostics/union_adaptive_20260923/CANONICAL_PLAN.md \
  --output "$UNION_RUN/canonical_preparation_replay"
```

These latter two commands are geometry/manifest operations, not scientific
execution. They select24 new sources and record4 exact completed pilot reuses.
All output paths must be new. The original exact submission wrappers/commands,
receipts and immutable snapshots remain inside canonical_proposals_v1 and
canonical_pool_v1. `prepare_pool --proposals <terminal proposal collection>
--agreement <same plan> --output <new directory>` prepares the existing shared
MACE/GFN2 runners; it does not execute them. No225-fold continuation or production
promotion is authorized by this canonical-only command list.
