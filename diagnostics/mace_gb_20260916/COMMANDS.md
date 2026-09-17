# Frozen-solvent runbook

Run from the repository root. Outputs are immutable; choose a new output
filename when collecting again. No new scientific inference is needed to
inspect the completed result.

```bash
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_hybrid.py dry-run --manifest workspaces/mace_gb_20260916/pilot_v2/manifest.json
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_hybrid.py collect --manifest workspaces/mace_gb_20260916/pilot_v2/manifest.json --output workspaces/mace_gb_20260916/recollection.json
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_hybrid.py report --collection workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --output workspaces/mace_gb_20260916/review.md
```

Prepare an exact scientific repeat in a new directory:

```bash
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_gb.py --medium-collection workspaces/mace_hydrogen_20260916/pilot_v2/medium/collection_job_1200681.json --large-collection workspaces/mace_hydrogen_20260916/pilot_v2/large/collection_job_1200682.json --core-collection workspaces/mace_analytic_20260916/pilot_v1/collection_job_1200470.json --agreement diagnostics/mace_gb_20260916/CUDA_REPAIR.md --output workspaces/mace_gb_20260916/repeat_prepared
```

The existing `run_pilot.sbatch` executes a finite manifest and reuses only
verified successful tasks. Its arguments are Python, absolute manifest path,
`native`, and optional explicit task IDs. For the A5000 node override partition
`gpu`, node `node-128-512g-8gpu-1`, `--gres=gpu:1`, `--cpus-per-task=16`,
`--mem=64474M`, and export `MACE_MIN_MEMORY_MIB=64474`. Failed attempts remain
visible. The scheduler's seven-day QOS limit is not a project compute budget.

The isolated environment uses OpenMM 8.5.1 plus CUDA12 8.5.1 wheels under
`workspaces/mace_gb_20260916/software_v2/`; download/install logs and wheels are
retained. Versions, native libraries and source hashes are in the prepared
software inventory. Do not install these over the shared Conda environment.
