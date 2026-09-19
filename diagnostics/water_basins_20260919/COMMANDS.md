# Finite rigid-water basin operations

Run from repository root. Completed jobs must not be rerun. Paths are explicit;
new preparation/collection/report outputs are write-once.

```bash
WATER_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
MACE_PY=$PWD/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

Read-only audit and tests:

```bash
"$WATER_PY" scripts/water_basin_sampling.py validate \
  --manifest workspaces/water_basins_20260919/sampling_v1/manifest.json
"$WATER_PY" -m unittest discover -s tests -p test_water_basin_sampling.py -v
"$MACE_PY" workspaces/water_basins_20260919/sampling_v1/implementation/mace_hybrid.py dry-run \
  --manifest workspaces/water_basins_20260919/sampling_v1/manifest.json
```

Sampling runs through the existing `mace_hybrid.py execute` and existing
`diagnostics/mace_hybrid_20260916/run_pilot.sbatch`, with native memory mode.
Its only additional dispatch is isolated inside the new manifest's implementation
snapshot. Shared scripts and previous frozen snapshots are unchanged.

Collection writes an explicit snapshot, including partial status:

```bash
"$WATER_PY" scripts/water_basin_sampling.py collect \
  --manifest workspaces/water_basins_20260919/sampling_v1/manifest.json \
  --output workspaces/water_basins_20260919/sampling_v1/collection_manual_v1.json
```

Each completed endpoint can prepare its fixed four native representatives without
waiting for the remaining endpoint sampling. Example for the first endpoint:

```bash
"$WATER_PY" scripts/water_basin_validation.py prepare \
  --collection workspaces/water_basins_20260919/sampling_v1/collection_manual_v1.json \
  --center-id 1F6S__11__Ca \
  --agreement diagnostics/water_basins_20260919/AGREEMENT.md \
  --output workspaces/water_basins_20260919/native_replay_1F6S_11_Ca_v1
"$WATER_PY" scripts/affordable_workflow.py dry-run \
  --manifest workspaces/water_basins_20260919/native_replay_1F6S_11_Ca_v1/manifest.json
sbatch --parsable \
  --output="$PWD/workspaces/water_basins_20260919/native_replay_1F6S_11_Ca_v1/slurm_%j.out" \
  --error="$PWD/workspaces/water_basins_20260919/native_replay_1F6S_11_Ca_v1/slurm_%j.err" \
  diagnostics/water_basins_20260919/run_dft.sbatch \
  "$PWD/workspaces/water_basins_20260919/native_replay_1F6S_11_Ca_v1/manifest.json"
```

The native job uses the existing allocation-aware ORCA manifest runner and retains
all actual output/analytic-gradient receipts. `water_basin_validation.py report`
accepts repeated `--validation` arguments for the four completed native collections
and one final `--sampling` collection. Actual completed collections are recorded below; no dummy successful execution
or occupancy values are supplied.


## Replay the actual completed comparison without compute

```bash
PYTHONPATH=workspaces/water_basins_20260919/sampling_v1/implementation \
"$WATER_PY" workspaces/water_basins_20260919/analysis_implementation_v1/water_basin_validation.py report \
  --sampling workspaces/water_basins_20260919/sampling_v1/collection_job_1202081.json \
  --validation workspaces/water_basins_20260919/native_1F6S_11_Ca_v2/collection_1202085.json \
  --validation workspaces/water_basins_20260919/native_1F6S_11_La_v1/collection_1202086.json \
  --validation workspaces/water_basins_20260919/native_6IP9_110_Ca_v1/collection_1202087.json \
  --validation workspaces/water_basins_20260919/native_6IP9_110_La_v1/collection_1202088.json \
  --output diagnostics/water_basins_20260919/export_replay_v1
```

Existing export_v1 is immutable. All 16 nativechecks completed; no scientific
retries occurred. The initial native_1F6S_11_Ca_v1 contains an unsubmitted
preparation interrupted by negligible weight-recomputation roundoff; preserve it.
