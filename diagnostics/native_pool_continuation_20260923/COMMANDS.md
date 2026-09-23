# Two-pass continuation operations

Run from the repository root with the pinned driver below. Actual submission:
job **1210333**, recorded in
`workspaces/native_pool_continuation_20260923/run_v1/SUBMISSION.json`.
The job performs both passes and collects even after an endpoint failure.
These commands describe existing artifacts; do not resubmit completed work.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/native_pool_continuation.py prepare \
  --inventory diagnostics/native_xtb_restart_20260922/SEED_AVAILABILITY.json \
  --agreement diagnostics/native_pool_continuation_20260923/PLAN.md \
  --output workspaces/native_pool_continuation_20260923/run_v1/stage1
```

The prepared manifest pins its executable implementation. The submitted wrapper
uses that implementation for execute/collect and prepares stage 2 from the exact
stage 1 collection. It then writes `run_v1/result.json`. If a job's final
reporting step is interrupted after collection, the current no-compute command
(including the corrected archived/continued receipt attribution) is:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/native_pool_continuation.py compare \
  --stage1 workspaces/native_pool_continuation_20260923/run_v1/stage1/collection.json \
  --stage2 workspaces/native_pool_continuation_20260923/run_v1/stage2/collection.json \
  --output workspaces/native_pool_continuation_20260923/run_v1/result_receipts_v2.json
```

This command requires a new output path if the result already exists. Neither
collection nor comparison launches calculations. Parser/actual-fixture tests:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  -m unittest discover -s tests -p test_native_pool_continuation.py -v
```

The existing ORCA runner is unchanged. `execute_tasks(manifest, validator)` is
the scoped seed-activation wrapper. Its optional `seeded=False` route is for the
separately declared derivative experiment's initial displaced calculations;
that caller must supply its own finite-scope validator and explicit NoAutostart.
The fixed-pool CLI always uses saved paired seeds. No extra scientific work is
authorized by the helper API.
