# Operations

The approved four new tasks are complete. These first commands only validate and
replay their real results. Run from the repository root; no source edits needed.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
WATER_REFERENCE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
WATER_REFERENCE_WORK="$PWD/workspaces/water_reference_validation_20260919"
"$WATER_REFERENCE_PY" scripts/water_reference_transfer.py validate \
  --manifest "$WATER_REFERENCE_WORK/proposals_v1/manifest.json"
"$WATER_REFERENCE_PY" -m unittest discover -s tests -p test_water_reference_transfer.py -v
```

Write-once collection replay (each output path must be absent):

```bash
"$WATER_REFERENCE_PY" "$WATER_REFERENCE_WORK/proposals_v1/implementation/water_reference_transfer.py" collect \
  --manifest "$WATER_REFERENCE_WORK/proposals_v1/manifest.json" \
  --output "$WATER_REFERENCE_WORK/proposals_v1/collection_review.json"
"$WATER_REFERENCE_PY" "$WATER_REFERENCE_WORK/dft_v1/implementation/water_reference_transfer.py" collect-dft \
  --manifest "$WATER_REFERENCE_WORK/dft_v1/manifest.json" \
  --output "$WATER_REFERENCE_WORK/dft_v1/collection_review.json"
```

## Preparation and dry-run interface

Explicit configuration supplies ordered sites, source-backed repaired cores,
baseline receipts, qualified preparation policy and agreement. The adapter can
prepare a new write-once reproduction manifest without launching calculations:

```bash
"$WATER_REFERENCE_PY" scripts/water_reference_transfer.py prepare \
  --config "$PWD/diagnostics/water_reference_validation_20260919/CONFIG.json" \
  --output "$WATER_REFERENCE_WORK/review_proposals_v1"
"$WATER_REFERENCE_PY" "$WATER_REFERENCE_WORK/review_proposals_v1/implementation/mace_hybrid.py" dry-run \
  --manifest "$WATER_REFERENCE_WORK/review_proposals_v1/manifest.json"
```

The original execution is recorded as exact argument arrays in
[MACE_SUBMISSION.json](MACE_SUBMISSION.json) and
[DFT_SUBMISSION.json](DFT_SUBMISSION.json); **do not resubmit completed tasks**.
MACE uses the existing `run_pilot.sbatch`, pinned native environment and manifest
closure. DFT uses [run_dft.sbatch](run_dft.sbatch), the existing bounded manifest
executor, two 16-rank endpoints and native ORCA. Accepted tasks remain cached and
partial failures remain visible.

The DFT bridge accepts an actual completed proposal collection and creates only
wet-site endpoint tasks; dry sites retain their original receipts. To inspect a
fresh preparation without rerunning scientific work:

```bash
"$WATER_REFERENCE_PY" scripts/water_reference_transfer.py prepare-dft \
  --collection "$WATER_REFERENCE_WORK/proposals_v1/collection_job_1202475.json" \
  --output "$WATER_REFERENCE_WORK/review_dft_v1"
"$WATER_REFERENCE_PY" scripts/affordable_workflow.py dry-run \
  --manifest "$WATER_REFERENCE_WORK/review_dft_v1/manifest.json"
```

Report: [REPORT.md](REPORT.md). The machine-readable paired table is
`workspaces/water_reference_validation_20260919/dft_v1/collection_1202476.json`.
This standalone experiment adapter does not replace the parent-owned production
contextual-water interface or enable unsupported whole-protein/cofactor inputs.
