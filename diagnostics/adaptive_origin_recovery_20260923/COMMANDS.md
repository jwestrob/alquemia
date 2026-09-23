# Direct-origin recovery commands

Run from the repository root. Jobs1210389/1210399 are complete; do not resubmit
their scientific work. Exact invocation receipts are in `run_v2/`.

Reproduce the no-compute result report with a fresh output path:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/adaptive_origin_recovery.py score \
  --collection workspaces/adaptive_origin_recovery_20260923/run_v2/pool/final_collection.json \
  --output workspaces/adaptive_origin_recovery_20260923/run_v2/decisions_replay.json
```

Actual source/method validation and real-fixture tests:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/adaptive_origin_recovery.py validate-pool \
  --manifest workspaces/adaptive_origin_recovery_20260923/run_v2/pool/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  -m unittest discover -s tests -p test_adaptive_origin_recovery.py -v
```

The scientific adapter exposes `prepare`, `execute`, `prepare-pool`,
`execute-cross`, `collect` and `score`, all with explicit paths and a fixed
two-source scope. `prepare-pool` reuses the existing native pool's low-level
manifest builder; `run_low.sbatch` invokes the existing ORCA executor and always
collects after execution failure. No new workflow or production entrypoint
is introduced. Full225 overlay belongs to the parent's minimal-pool report.
