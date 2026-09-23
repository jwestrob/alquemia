# Twenty-call strict native stopping diagnostic

Job1211229 executes the finite plan. Do not duplicate its calculations or add
an automatic restart. SUBMISSION.json pins the exact command, wrapper, manifest
and agreement. It always attempts collection after execution errors.

Read-only validation:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/strict_native_stopping_20260923/run_v1/implementation/strict_native_stopping.py \
validate --manifest workspaces/strict_native_stopping_20260923/run_v1/manifest.json
```

Replay collection to a new output after actual execution:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/strict_native_stopping_20260923/run_v1/implementation/strict_native_stopping.py \
collect --manifest workspaces/strict_native_stopping_20260923/run_v1/manifest.json \
--output workspaces/strict_native_stopping_20260923/run_v1/COLLECTION_replay_v1.json
```

No full classifier score or reference is defined from these ten incomplete
metal/medium comparisons. Every source, output, seed and failure remains
available, and neither initialization is selected by its resulting energy.
