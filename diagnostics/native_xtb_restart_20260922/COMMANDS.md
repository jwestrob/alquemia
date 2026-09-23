# Exact eight-test restart operation and activation recovery

Executed preparation is `workspaces/native_xtb_restart_20260922/run_v2/manifest.json`,
SHA `951014e426902073f87fbd5b7d0eeff337b4c5e43e5e3b4e0ea43befd47914d2`.
Job1210179 completed all eight calls and automatic collection. Do not submit
it again. `run_v1` is a preparation-only predecessor that was never submitted;
v2 adds the shared explicit user-authorization pin without changing scientific settings.

The approved activation recovery also completed, job1210185, with all eight native
restarts confirmed. Its manifest is `recovery_v1/manifest.json` under the same
workspace, SHA `756017eafc7d1762ae0869081aad7882aafaf2c29d34a4901e7d640b3999ab38`.
Sixteen actual calls in total include the original unsuccessful activation.

From the repository root, use the pinned driver:

```bash
RESTART_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$RESTART_PY" -m unittest discover -s tests -p test_native_xtb_restart.py -v
```

The immutable original snapshot is
`workspaces/native_xtb_restart_20260922/run_v2/collection_1210179.json`.
The recovery snapshot is `recovery_v1/collection_1210185.json`.
For a read-only recovery collection replay to a fresh output path:

```bash
"$RESTART_PY" workspaces/native_xtb_restart_20260922/recovery_v1/implementation/native_xtb_restart.py collect \
  --manifest workspaces/native_xtb_restart_20260922/recovery_v1/manifest.json \
  --output workspaces/native_xtb_restart_20260922/recovery_v1/collection_replay.json
```

This starts no calculation. Each restart requires an immutable seed hash, actual
`XTBRESTART` output, native method/state, parameters, charge checks and the existing
executor's complete receipt. Unconfirmed/failed cells retain explicit statuses;
no fallback replaces them. The raw observed energy can be retained even when a
restart cannot be confirmed; it is distinct from the qualified matrix energy.

All exact launch arguments and wrapper pins are in each run's `SUBMISSION.json`.
`SEEDS_BEFORE.json`, `SEEDS_AFTER.json`, per-task immutable/active/after seed files,
outputs and execution receipts are retained. No scratch cleanup touches these.
The recovery additionally retains each matching archived GBW, its active/post-run
copy and AutoStart GES. See `REPORT.md` for the result and qualification limits.
