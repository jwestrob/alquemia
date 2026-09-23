# Exact bounded recovery commands

Run from the repository root. Job1211825 already executes the two declared calls;
do not resubmit completed attempts. Production and the primary full100 result
stay unchanged.

```bash
PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RUN="$PWD/workspaces/native_failed_cell_recovery_20260923/run_v1"
"$PY" "$RUN/implementation/native_failed_cell_recovery.py" validate \
  --manifest "$RUN/manifest.json"
"$PY" -m unittest discover -s tests -p test_native_failed_cell_recovery.py -v
```

The immutable submitted wrapper executes and always collects, preserving failed
receipts. Collection is write-once at `run_v1/COLLECTION.json`; no molecular calls
occur in validation, collection or reporting. Exact `SUBMISSION.json` retains
command, manifest and wrapper pins. `SEEDS_BEFORE/AFTER.json` retains both seed
identities at the actual runtime basename. `EXECUTION.json` records solver wall
time and job identity, separate from full allocation accounting.

`status=qualified_recovery_sensitivity` requires both real native restarts and
≤0.1 kcal/mol agreement. Its sole proposed energy is the origin-seeded energy;
it is never selected by lower energy or class. Otherwise energy remains null.
`primary_status=unavailable_unchanged` applies regardless of the diagnostic result.
