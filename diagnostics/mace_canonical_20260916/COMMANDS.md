# Canonical MACE candidate operations

Run from the repository root. PLAN.md is frozen before new results. Every output
writer refuses overwrite. Read AUDIT.md for legacy serialization and grouping.
These commands preserve the production baseline and launch no automatic rescore.

## Audit a fresh inventory

```bash
python scripts/mace_canonical.py audit --calibration diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json --holdouts diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/result/holdout_result.json --external-preparation diagnostics/pqq_q46444_1kb0_external_validation_20260915/prepared/external_preparation.json --external-collection workspaces/baseline_benchmark_20260915/run_v1/collection_1199299.json --panel testset_expansion/frozen_v0_ids.tsv --agreement diagnostics/mace_canonical_20260916/PLAN.md --output workspaces/mace_canonical_20260916/audit_replay_v1
```

## Prepare and validate exact MACE tasks

```bash
python scripts/mace_canonical_run.py prepare-mace --source-inventory workspaces/mace_canonical_20260916/audit_v2/inventory.json --parent workspaces/mace_global_benchmark_20260916/mace_v1/medium/collection_job_1200701.json --checkpoint medium --reference-mace workspaces/mace_local_correction_20260916/mace_v2/collection_job_1200717.json --reference-gb workspaces/mace_local_correction_20260916/gb_v1/collection_job_1200719.json --agreement diagnostics/mace_canonical_20260916/PLAN.md --output workspaces/mace_canonical_20260916/mace_replay_v1/medium
python scripts/mace_canonical_run.py prepare-mace --source-inventory workspaces/mace_canonical_20260916/audit_v2/inventory.json --parent workspaces/mace_global_benchmark_20260916/mace_v1/large/collection_job_1200702.json --checkpoint large --agreement diagnostics/mace_canonical_20260916/PLAN.md --output workspaces/mace_canonical_20260916/mace_replay_v1/large
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_canonical_20260916/mace_v2/medium/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_canonical_20260916/mace_v2/medium/manifest.json
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_canonical_20260916/mace_v2/large/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_canonical_20260916/mace_v2/large/manifest.json
```

## Execute / recover

Actual submitted argv, IDs and manifest hashes are in each mace_v2 checkpoint's
submission.json:1200776medium and1200779large. OneA5000,16CPU,64474MiB per job.
Existing run_pilot.sbatch executes the manifest's copied implementation and
collects on exit. The seven-day scheduler QOS request is not a project spending
budget. Do not submit concurrent executors on one manifest or remove locks.
Valid cached attempts are reused; failed/partial attempts remain visible.
For the fresh medium replay above, the complete execution command is:

```bash
ALQUEMIA_CANONICAL_ROOT=$(pwd -P)
sbatch --parsable --partition=gpu --nodelist=node-128-512g-8gpu-1 --job-name=alquemia_canonical_replay --cpus-per-task=16 --gres=gpu:1 --mem=64474M --time=7-00:00:00 --export=ALL,MACE_MIN_MEMORY_MIB=64474 --output="$ALQUEMIA_CANONICAL_ROOT/workspaces/mace_canonical_20260916/mace_replay_v1/medium/slurm_%j.out" --error="$ALQUEMIA_CANONICAL_ROOT/workspaces/mace_canonical_20260916/mace_replay_v1/medium/slurm_%j.err" diagnostics/mace_hybrid_20260916/run_pilot.sbatch "$ALQUEMIA_CANONICAL_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" "$ALQUEMIA_CANONICAL_ROOT/workspaces/mace_canonical_20260916/mace_replay_v1/medium/manifest.json" native
```

## Prepare solvent after actual MACE completion

```bash
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_canonical_run.py prepare-gb --collection workspaces/mace_canonical_20260916/mace_v2/medium/collection_job_1200776.json --solver-validation workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --agreement diagnostics/mace_canonical_20260916/PLAN.md --output workspaces/mace_canonical_20260916/gb_replay_v1/medium
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_canonical_run.py prepare-gb --collection workspaces/mace_canonical_20260916/mace_v2/large/collection_job_1200779.json --solver-validation workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --agreement diagnostics/mace_canonical_20260916/PLAN.md --output workspaces/mace_canonical_20260916/gb_replay_v1/large
```

GB execution uses the same wrapper and resource allocation, with the literal
absolute workspaces/mace_gb_20260916/software_v2/venv/bin/python path. Do not
resolve the venv symlink or reuse the MACE Python for CUDA12OpenMM. Actual GB
submission/collection paths are alongside their finite manifests: jobs 1200781
(medium) and 1200782 (large), both complete.

## Collect / compare / report

Both actual gb_v1 jobs finished. These read-only collectors produce explicit
report inputs; no scientific inference is launched:

```bash
python scripts/mace_hybrid.py collect --manifest workspaces/mace_canonical_20260916/gb_v1/medium/manifest.json --output workspaces/mace_canonical_20260916/medium_gb_report_collection.json
python scripts/mace_hybrid.py collect --manifest workspaces/mace_canonical_20260916/gb_v1/large/manifest.json --output workspaces/mace_canonical_20260916/large_gb_report_collection.json
python scripts/mace_canonical_report.py --medium-mace workspaces/mace_canonical_20260916/mace_v2/medium/collection_job_1200776.json --medium-gb workspaces/mace_canonical_20260916/medium_gb_report_collection.json --large-mace workspaces/mace_canonical_20260916/mace_v2/large/collection_job_1200779.json --large-gb workspaces/mace_canonical_20260916/large_gb_report_collection.json --output workspaces/mace_canonical_20260916/report_replay_v1
python -m unittest discover -s tests -p test_mace_canonical.py -v
```

The report verifies receipts, calibrates only designated rows, retains unavailable
cases and all denominators, and compares checkpoints without selecting a winner.
No aquo reference or old threshold is inherited. Raw component contrasts and
published baseline results stay side by side; whole-family transfer limitations
and functional-versus-affinity distinctions remain explicit.
