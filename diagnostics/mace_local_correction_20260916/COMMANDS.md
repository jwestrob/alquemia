# Matched local correction runbook

Run from the repository root. Scope is PLAN.md. Old XYZs and all ten archived
DFT endpoints were verified against successful execution receipts. New global-H
geometry has no DFT endpoint yet and must not be reported as a valid hybrid.

Prepared executable manifest: `workspaces/mace_local_correction_20260916/mace_v2/manifest.json`.
Job1200717 executes20 medium MACE calls. The unexecuted v1 preserves initial
preparation; v2 has byte-identical scientific inputs and adds tolerance for
last-bit cross-CPU mapping-norm roundoff plus stricter source-receipt checks.
Four real-fixture preparation/parser/recovery tests pass.

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python scripts/mace_hybrid.py dry-run --manifest workspaces/mace_local_correction_20260916/mace_v2/manifest.json
```

After the MACE collection is complete, prepare the declared20-call solvent stage:

```bash
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_local_correction.py prepare-gb --collection workspaces/mace_local_correction_20260916/mace_v2/collection_job_1200717.json --global-gb workspaces/mace_global_benchmark_20260916/gb_v1/medium/collection_job_1200711.json --solver-validation workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --agreement diagnostics/mace_local_correction_20260916/PLAN.md --output workspaces/mace_local_correction_20260916/gb_v1
```

Use the existing runner with the pinned CUDA12 Python. This runnable command
submits the finite solvent manifest after preparation:

```bash
ALQUEMIA_LOCAL_ROOT=$(pwd -P)
sbatch --parsable --partition=gpu --nodelist=node-128-512g-8gpu-1 --job-name=alquemia_local_GB --cpus-per-task=16 --gres=gpu:1 --mem=64474M --time=7-00:00:00 --export=ALL,MACE_MIN_MEMORY_MIB=64474 --output="$ALQUEMIA_LOCAL_ROOT/workspaces/mace_local_correction_20260916/gb_v1/slurm_%j.out" --error="$ALQUEMIA_LOCAL_ROOT/workspaces/mace_local_correction_20260916/gb_v1/slurm_%j.err" diagnostics/mace_hybrid_20260916/run_pilot.sbatch "$ALQUEMIA_LOCAL_ROOT/workspaces/mace_gb_20260916/software_v2/venv/bin/python" "$ALQUEMIA_LOCAL_ROOT/workspaces/mace_local_correction_20260916/gb_v1/manifest.json" native
```

The EXIT handler collects success/failure. To collect independently, specify a
fresh output path; no failed correction is replaced with zero:

```bash
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_hybrid.py collect --manifest workspaces/mace_local_correction_20260916/gb_v1/manifest.json --output workspaces/mace_local_correction_20260916/solvent_review.json
```
