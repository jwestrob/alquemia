# Replay operations

Run from the repository root. Job1211311 is complete; do not resubmit.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RUN="$PWD/workspaces/strict_native_comparator_20260923/run_v1"
"$PY" "$RUN/implementation/strict_native_comparator.py" validate --manifest "$RUN/manifest.json"
"$PY" -m unittest discover -s tests -p test_strict_native_comparator.py -v
```

The existing immutable `COLLECTION.json` contains `rows` (all48 scalar receipts) and `cases` (all4 matrices, complete pools and frozen-reference calls). Read these directly for matched comparison; no source edits or new calculation is needed. `COSTS.json` retains allocation and engine runtime.
