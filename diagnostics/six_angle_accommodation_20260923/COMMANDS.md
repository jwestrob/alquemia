# Six-angle pilot operations

All paths are explicit. `run_v2` is the executed immutable preparation;
`run_v1` failed preflight on the older origin-receipt schema and was never run.
The technical fix uses the existing archive parser; no source or science changed.

## Prepared and submitted

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
scripts/six_angle_accommodation.py prepare \
  --feasibility workspaces/c5ax_response_20260923/SIX_MODE_FEASIBILITY_v1.json \
  --reference workspaces/strict_native_pool_20260923/run_v1/REFERENCES.json \
  --agreement diagnostics/six_angle_accommodation_20260923/PLAN.md \
  --output workspaces/six_angle_accommodation_20260923/run_v2
```

The output is write-once; do not rerun this command at an existing path.
Actual manifest SHA256 is
`7c59b163736462a1d3ffc7643b5b667cb67531e220ff1879c8ee3e0bf8d87d8a`.
GPU job1211537 uses `run_gpu.sbatch`. Its molecular work completed; subsequent
scalar staging failed because the audit adapter needed the complete archived
task metadata. Both empty staging attempts remain. Dependency1211538 never ran
and was canceled after confirming PENDING/zero runtime. Corrected scalar job
1211603 uses the same `run_scalar.sbatch` and `run_v2/scalar_v3/manifest.json`,
SHA `f351819c6703c245e10258bf3674cf91489ee0b0efcdbff06e0979e3eeb4aed6`.
The JSON submission/replacement receipts retain the exact resources and
scheduler records. No MACE calculation was repeated; do not duplicate the jobs.

## Read-only preflight

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/six_angle_accommodation_20260923/run_v2/implementation/six_angle_accommodation.py \
  validate --manifest workspaces/six_angle_accommodation_20260923/run_v2/manifest.json
```

The prepared checks replay all18warm starts and paired-origin selections. The
five real-fixture tests in `TESTS_preflight_v2.txt` also check six-column physical
constraint derivatives, Cartesian-gradient work and strict missing-cell handling.
These tests use archived scientific results, not new model evaluations.

## Actual stages and collection

The GPU wrapper calls `execute`, then `prepare_low`; it writes
`run_v2/GPU_COLLECTION.json`; repaired staging writes
`run_v2/scalar_v3/manifest.json`. This stays within
the finite18searches/18cross/72scalar ceiling. Failed proposals remain in the
source collection; no automatic alternate starts are launched.

The CPU wrapper calls `execute_low`, then `collect`, writing
`run_v2/scalar_v3/{EXECUTION,COLLECTION}.json`. A read-only recollection, if needed
after completed execution, must use a fresh output filename:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/six_angle_accommodation_20260923/run_v2/scalar_v3/implementation/six_angle_accommodation.py \
  collect --manifest workspaces/six_angle_accommodation_20260923/run_v2/scalar_v3/manifest.json \
  --output workspaces/six_angle_accommodation_20260923/run_v2/scalar_v3/COLLECTION_replay.json
```

Interpret per-source `own_added_work`, `old_band_transfer`, full five-geometry
`matrix` and mathematical/operational `pool`. These are development results
under frozen strict32 bands. No reference is fitted from the nine sources.
