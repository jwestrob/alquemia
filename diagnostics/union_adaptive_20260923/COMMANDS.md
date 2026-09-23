# Consistent-context × adaptive research pilot

The declared scope is eight consumed sources, sixteen physical searches and the
same three-candidate pool for both metals. No production/default change or new
reference. Existing completed molecular jobs should not be resubmitted.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
UNION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
UNION_RUN=$PWD/workspaces/union_adaptive_20260923
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$UNION_PY" scripts/union_adaptive.py validate --manifest "$UNION_RUN/proposals_v1/manifest.json"
"$UNION_PY" scripts/nikasha_pool.py validate --manifest "$UNION_RUN/pool_v1/manifest.json"
"$UNION_PY" -m unittest discover -s tests -p test_union_adaptive.py -v
```

Read-only comparison after the finite jobs and dependent collector finish:

```bash
"$UNION_PY" scripts/union_adaptive.py compare \
  --collection "$UNION_RUN/pool_v1/collection_final.json" \
  --union_reference workspaces/consistent_context_20260922/calibration28_v1/REFERENCE_v2.json \
  --output "$UNION_RUN/comparison_replay.json"
```

The output path must be new. Exact execution pins and wrapper commands are in
proposals_v1/SUBMISSION.json, pool_v1/SUBMISSIONS.json and
pool_v1/SUBMISSION_PLACEMENT_v2.json. Pending-only solver/collector requests on
node128 were replaced on node224 without changing scientific input; both sets
remain recorded. The latter collector depends on both completed cross-MACE and
final GFN2 jobs and writes `pool_v1/collection_final.json`, retaining failures.

Geometry-only preparation is fully parameterized and writes a new workspace:

```bash
"$UNION_PY" scripts/union_adaptive.py prepare \
  --preparation workspaces/consistent_context_20260922/prepared_v1/PREPARATION.json \
  --crystals workspaces/consistent_context_20260922/crystal_controls_v1/CRYSTALS.json \
  --calibration workspaces/consistent_context_20260922/calibration28_v1/collection_final_v2.json \
  --transfer workspaces/consistent_context_20260922/primary225_v1/collection_final_v1.json \
  --inputs diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json \
  --source workspaces/adaptive_completion_20260922/original30_v1/manifest.json \
  --agreement diagnostics/union_adaptive_20260923/PLAN.md \
  --output "$UNION_RUN/preparation_replay"
```

`execute --manifest` requires the existing32CPU/200000MiB GPU allocation;
`collect --manifest --output` only parses actual proposal receipts. `prepare_pool`
accepts `--proposals`, `--agreement`, `--output`, then uses unchanged shared
`nikasha_pool.py execute-mace/collect` and existing native ORCA executors. The new
shared dispatch validates only this explicitly named research protocol.
