# Compact matched solvent transfer: recorded operations

These commands use the approved, immutable 33-case experiment. They do not change
the released DFT/water-preparation baseline. Source records, method and selection
are in `PLAN.md`, `INVENTORY.json` and the three submission receipts.

```bash
COMPACT_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
COMPACT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$COMPACT_ROOT"
```

## Audit and tests: no scientific execution

```bash
"$COMPACT_PY" scripts/compact_solvation.py validate \
  --manifest "$COMPACT_ROOT/workspaces/compact_solvation_20260920/full_v1/manifest.json"
"$COMPACT_PY" -m unittest discover -s tests -p 'test_compact_solvation*.py' -v
```

Validation verifies the actual source, code, executable, charge/spin, coordinates,
input recipe, finite task inventory and reused calculation receipts. It includes
the existing runner's dry-run validation and launches no ORCA calculation.

## Actual preparation and execution already performed

The following are the exact full-panel preparation/submission operations. Their
outputs already exist; **do not resubmit completed calculations**. Preparation
refuses an existing destination. The submission receipt records job1203165.

```bash
"$COMPACT_PY" scripts/compact_solvation.py prepare \
  --inventory "$COMPACT_ROOT/diagnostics/compact_solvation_20260920/INVENTORY.json" \
  --agreement "$COMPACT_ROOT/diagnostics/compact_solvation_20260920/PLAN.md" \
  --output "$COMPACT_ROOT/workspaces/compact_solvation_20260920/full_v1" \
  --reuse-manifest "$COMPACT_ROOT/workspaces/compact_solvation_20260920/pilot_v1/manifest.json"
sbatch --parsable \
  --output="$COMPACT_ROOT/workspaces/compact_solvation_20260920/full_v1/slurm_%j.out" \
  --error="$COMPACT_ROOT/workspaces/compact_solvation_20260920/full_v1/slurm_%j.err" \
  "$COMPACT_ROOT/diagnostics/compact_solvation_20260920/run.sbatch" \
  "$COMPACT_ROOT/workspaces/compact_solvation_20260920/full_v1/manifest.json"
```

The batch uses 64CPUs, eight concurrent eight-rank tasks, 128GiB and no GPU.
Its private implementation snapshot executes exactly256 new GFN2 tasks and reuses
eight pilot results. There is no additional DFT/MACE call, optimization or project
time budget. Failed attempts remain visible. The fixed eight alternate-SCF
checks failed and are not retried by this manifest.

## Collect and compare actual output

Collection and comparison do no scientific calculation. Existing output paths
are immutable: these commands refuse to overwrite them. The batch collects after
successful execution; manual collection remains available for partial failure.

```bash
"$COMPACT_PY" scripts/compact_solvation.py collect \
  --manifest "$COMPACT_ROOT/workspaces/compact_solvation_20260920/full_v1/manifest.json" \
  --output "$COMPACT_ROOT/workspaces/compact_solvation_20260920/full_v1/collection_1203165.json"
"$COMPACT_PY" scripts/compact_solvation_compare.py compare \
  --inventory "$COMPACT_ROOT/diagnostics/compact_solvation_20260920/INVENTORY.json" \
  --collection "$COMPACT_ROOT/workspaces/compact_solvation_20260920/full_v1/collection_1203165.json" \
  --numerical-collection "$COMPACT_ROOT/workspaces/compact_solvation_20260920/numerical_v1/failed_collection_v1.json" \
  --solver native \
  --output "$COMPACT_ROOT/workspaces/compact_solvation_20260920/full_v1/comparison_v1.json"
```

Missing/failed endpoints remain unavailable. The new comparison has its own
calibration; the failed numerical check is retained even if discrimination helps.
An aquo-referenced score is unavailable. Nothing substitutes old scores for failed
new calculations.

## Read results and actual allocation cost

```bash
cat diagnostics/compact_solvation_20260920/REPORT.md
sacct -n -P -j 1203145,1203147,1203165 \
  --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,TotalCPU,MaxRSS,ReqMem,AllocTRES,NodeList
```

The report separates added GFN2 cost from reused native MACE inference and from
one-time numerical checks. No matched end-to-end scanner speedup is inferred from
cached MACE results or mixed GPU hardware.
