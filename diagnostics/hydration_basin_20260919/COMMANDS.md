# Coupled physical water response

Run from the repository root. These operations use the existing MACE and ORCA
runners. Preparation and reports are write-once; keep failed attempts. The default
scorer is unchanged. Do not resubmit completed scientific calculations.

```bash
HYDRATION_ROOT=$PWD
HYDRATION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
MACE_PY="$HYDRATION_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Audit existing manifests and actual results

```bash
"$HYDRATION_PY" scripts/hydration_basin.py collect \
  --manifest workspaces/hydration_basin_20260919/mace_v1/manifest.json \
  --output workspaces/hydration_basin_20260919/mace_v1/recollection_v1.json
"$HYDRATION_PY" scripts/affordable_workflow.py dry-run \
  --manifest workspaces/hydration_basin_20260919/validation_v1/dft/manifest.json
"$MACE_PY" workspaces/hydration_basin_20260919/validation_v1/mace/implementation/mace_hybrid.py dry-run \
  --manifest workspaces/hydration_basin_20260919/validation_v1/mace/manifest.json
```

## Prepare a fresh replay, without submitting it

The first operation prepares all 20 nonempty states from the actual completed
occupancy table. Native validation needs the resulting *completed* MACE collection.
Do not substitute a manifest or predicted energy for that collection.

```bash
"$HYDRATION_PY" scripts/hydration_basin.py prepare \
  --collection workspaces/hydration_occupancy_20260918/dft_v1/collection_1201910.json \
  --inventory workspaces/mace_canonical_20260916/audit_v2/inventory.json \
  --agreement diagnostics/hydration_basin_20260919/AGREEMENT.md \
  --output workspaces/hydration_basin_20260919/replay_mace_v1
"$MACE_PY" workspaces/hydration_basin_20260919/replay_mace_v1/implementation/mace_hybrid.py dry-run \
  --manifest workspaces/hydration_basin_20260919/replay_mace_v1/manifest.json
```

This separate operation prepares a replay of the actual 20 native validation
geometries from the completed Stage A calculation:

```bash
"$HYDRATION_PY" scripts/hydration_basin_validation.py prepare \
  --collection workspaces/hydration_basin_20260919/mace_v1/collection_job_1202050.json \
  --agreement diagnostics/hydration_basin_20260919/BOUNDED_CONTINUATION.md \
  --output workspaces/hydration_basin_20260919/replay_validation_v1
```

## Execute a prepared finite manifest

These commands target the fresh replay paths above. They are documented operations,
not a request to repeat the completed experiment. Use an available authorized GPU
allocation; preserve the recorded memory fraction and measure actual resources.

```bash
sbatch --parsable --partition=gpu --nodelist=node-224-2t-8gpu-1 \
  --job-name=water-coupled-mace --cpus-per-task=16 --gres=gpu:1 --mem=200000M \
  --export=ALL,MACE_MIN_MEMORY_MIB=200000 \
  --output="$HYDRATION_ROOT/workspaces/hydration_basin_20260919/replay_mace_v1/slurm_%j.out" \
  --error="$HYDRATION_ROOT/workspaces/hydration_basin_20260919/replay_mace_v1/slurm_%j.err" \
  diagnostics/mace_hybrid_20260916/run_pilot.sbatch "$MACE_PY" \
  "$HYDRATION_ROOT/workspaces/hydration_basin_20260919/replay_mace_v1/manifest.json" native
sbatch --parsable \
  --output="$HYDRATION_ROOT/workspaces/hydration_basin_20260919/replay_validation_v1/dft/slurm_%j.out" \
  --error="$HYDRATION_ROOT/workspaces/hydration_basin_20260919/replay_validation_v1/dft/slurm_%j.err" \
  diagnostics/hydration_basin_20260919/run_dft.sbatch \
  "$HYDRATION_ROOT/workspaces/hydration_basin_20260919/replay_validation_v1/dft/manifest.json"
sbatch --parsable --partition=gpu --nodelist=node-224-2t-8gpu-1 \
  --job-name=water-check-mace --cpus-per-task=16 --gres=gpu:1 --mem=200000M \
  --export=ALL,MACE_MIN_MEMORY_MIB=200000 \
  --output="$HYDRATION_ROOT/workspaces/hydration_basin_20260919/replay_validation_v1/mace/slurm_%j.out" \
  --error="$HYDRATION_ROOT/workspaces/hydration_basin_20260919/replay_validation_v1/mace/slurm_%j.err" \
  diagnostics/mace_hybrid_20260916/run_pilot.sbatch "$MACE_PY" \
  "$HYDRATION_ROOT/workspaces/hydration_basin_20260919/replay_validation_v1/mace/manifest.json" native
```

## Recenter and recheck

Recentering requires actual completed native energy/gradient checks. The original response
direction and each actual proposal-energy check must pass. Failed soft-mode checks
remain entropy failures. The next manifest includes only explicitly eligible states.
No numerical DFT gradients or unchecked cheap correction substitutes for missing
calculations. The finite continuation allows two rounds, keeping the same physical
trust region; it is not an open-ended native geometry optimization.

`hydration_basin_validation.py prepare --proposals-only` prepares actual native
checks for a completed recentered calculation. It does not pretend to repeat the
original directional curvature validation. Missing direction checks remain null.

These concrete preparation-only replays preserve the actual agreed selections:

```bash
"$HYDRATION_PY" scripts/hydration_basin_recenter.py \
  --validation diagnostics/hydration_basin_20260919/VALIDATION_v1.json \
  --agreement diagnostics/hydration_basin_20260919/BOUNDED_CONTINUATION.md \
  --output workspaces/hydration_basin_20260919/recenter_replay_g1
"$HYDRATION_PY" scripts/hydration_basin_validation.py prepare --proposals-only \
  --collection workspaces/hydration_basin_20260919/recenter_g1/collection_job_1202056.json \
  --agreement diagnostics/hydration_basin_20260919/BOUNDED_CONTINUATION.md \
  --output workspaces/hydration_basin_20260919/recenter_replay_g1_validation
```

## Compare and export the actual completed experiment

The initial native collection was recovered by parsing existing outputs; no DFT
rerun occurred. See [collection recovery](COLLECTION_RECOVERY.md). The frozen
analyzer reproduces all native comparisons. All paths below refer to actual jobs.

```bash
HYDRATION_ANALYZER=workspaces/hydration_basin_20260919/analysis_implementation_v4
"$HYDRATION_PY" "$HYDRATION_ANALYZER/hydration_basin_validation.py" compare \
  --preparation workspaces/hydration_basin_20260919/validation_v1/preparation.json \
  --mace-collection workspaces/hydration_basin_20260919/validation_v1/mace/collection_job_1202053.json \
  --dft-collection workspaces/hydration_basin_20260919/validation_v1/dft/collection_recovered_v1.json \
  --output diagnostics/hydration_basin_20260919/VALIDATION_replay_v1.json
"$HYDRATION_PY" "$HYDRATION_ANALYZER/hydration_basin_validation.py" compare \
  --preparation workspaces/hydration_basin_20260919/recenter_g1_validation/preparation.json \
  --mace-collection workspaces/hydration_basin_20260919/recenter_g1_validation/mace/collection_job_1202059.json \
  --dft-collection workspaces/hydration_basin_20260919/recenter_g1_validation/dft/collection_1202058.json \
  --output diagnostics/hydration_basin_20260919/VALIDATION_replay_g1.json
"$HYDRATION_PY" "$HYDRATION_ANALYZER/hydration_basin_validation.py" compare \
  --preparation workspaces/hydration_basin_20260919/recenter_g2_validation/preparation.json \
  --mace-collection workspaces/hydration_basin_20260919/recenter_g2_validation/mace/collection_job_1202065.json \
  --dft-collection workspaces/hydration_basin_20260919/recenter_g2_validation/dft/collection_1202064.json \
  --output diagnostics/hydration_basin_20260919/VALIDATION_replay_g2.json
"$HYDRATION_PY" "$HYDRATION_ANALYZER/hydration_basin_report.py" \
  --stage-a workspaces/hydration_basin_20260919/mace_v1/collection_job_1202050.json \
  --validation diagnostics/hydration_basin_20260919/VALIDATION_v1.json \
  --validation diagnostics/hydration_basin_20260919/VALIDATION_g1.json \
  --validation diagnostics/hydration_basin_20260919/VALIDATION_g2.json \
  --accounting diagnostics/hydration_basin_20260919/THERMODYNAMIC_ACCOUNTING_v2.md \
  --water-reference workspaces/hydration_network_20260918/water_reference_v2/reference_1201831.json \
  --output diagnostics/hydration_basin_20260919/export_replay_v1
```

Reports retain all failed checks and show actual Ca−La electronic contributions,
latest native geometries with explicit endpoint reuse, physical geometry checks,
and conditional thermal-extent diagnostics. Missing entropy/occupancy stays null.
Round2 requested48CPUs/192GiB for3concurrent16-rank tasks; Slurm allocated the
64CPU node, which is what the cost record charges. The already stationary6IP9La
endpoint was reused from round1. No extra chemistry was run to fill idle cores.

## Regression tests

```bash
"$HYDRATION_PY" -m unittest discover -s tests -p 'test_hydration*.py' -v
"$HYDRATION_PY" -m unittest discover -s tests -p test_affordable_development.py -v
```

Coordinate algebra tests use real atom sets and gradients. Native integration
results come only from the actual pinned MACE/ORCA receipts. No synthetic energies,
dummy structures, invented curvature or occupancy labels are used.
