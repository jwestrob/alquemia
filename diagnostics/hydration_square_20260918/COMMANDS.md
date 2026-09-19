# Hydration-square operations

Run from the repository root below. Completed inputs/outputs are immutable;
choose a new output name for an additional replay rather than overwriting it.
These commands use real pinned local artifacts. Production is unchanged.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
HYDRATION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
HYDRATION_RUN=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/hydration_square_20260918/repaired_v1
```

## Inspect and replay without electronic calculations

```bash
"$HYDRATION_PY" scripts/hydration_square.py audit \
  --config diagnostics/hydration_square_20260918/EXECUTION_CONFIG.json \
  --output workspaces/hydration_square_20260918/source_audit_replay1.json
"$HYDRATION_PY" scripts/affordable_workflow.py dry-run \
  --manifest "$HYDRATION_RUN/manifest.json" \
  --output "$HYDRATION_RUN/dry_run_replay1.json"
PYTHONPATH="$HYDRATION_RUN/implementation_1201801" "$HYDRATION_PY" \
  "$HYDRATION_RUN/implementation.py" collect \
  --manifest "$HYDRATION_RUN/manifest.json" \
  --output "$HYDRATION_RUN/collection_replay1.json"
```

Collection reports missing/failed endpoints explicitly if invoked before all
calculations finish. It never substitutes an archived full state for a repaired
state. Component SCF includes CPCM; do not add SCF and CPCM together.

## Recorded preparation and submission

These exact operations already prepared and submitted job1201801. Do not
repeat submission while it is running or after successful completion.

```bash
"$HYDRATION_PY" scripts/hydration_square.py prepare \
  --config diagnostics/hydration_square_20260918/EXECUTION_CONFIG.json \
  --output "$HYDRATION_RUN"
sbatch --parsable \
  --output="$HYDRATION_RUN/slurm_%j.out" \
  --error="$HYDRATION_RUN/slurm_%j.err" \
  diagnostics/hydration_square_20260918/run.sbatch "$HYDRATION_RUN/manifest.json"
```

The runner checks immutable input hashes, snapshots its implementation and uses
four concurrent 16-rank endpoints in a 64-CPU allocation. All14 endpoints belong
to the finite scientific manifest; no artificial compute-budget stop is imposed.
Failed/partial attempts are preserved and require an explicit fresh retry path.
No water occupancy, bulk-water chemical potential or calibrated classification
is produced by this diagnostic.
