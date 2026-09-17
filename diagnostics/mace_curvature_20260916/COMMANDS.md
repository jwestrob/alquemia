# Curvature pilot operations

Run from the repository root. PLAN.md fixes the scientific inputs and gates.
Preparation reuses exact archived DFT geometries and successful receipts; no DFT
is launched. Every operation refuses to overwrite an existing output.

## Audit / prepare / dry-run

The executed inputs are `workspaces/mace_curvature_20260916/mace_v1/medium/manifest.json`
and the corresponding `large/manifest.json`. Rebuild into a fresh named directory:

```bash
python scripts/mace_curvature.py prepare-mace --dft workspaces/ggr_mechanism_20260915/report_v1/collection_c.json --parent workspaces/mace_global_benchmark_20260916/mace_v1/medium/collection_job_1200701.json --agreement diagnostics/mace_curvature_20260916/PLAN.md --output workspaces/mace_curvature_20260916/replay/medium
python scripts/mace_curvature.py prepare-mace --dft workspaces/ggr_mechanism_20260915/report_v1/collection_c.json --parent workspaces/mace_global_benchmark_20260916/mace_v1/large/collection_job_1200702.json --agreement diagnostics/mace_curvature_20260916/PLAN.md --output workspaces/mace_curvature_20260916/replay/large
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python scripts/mace_hybrid.py dry-run --manifest workspaces/mace_curvature_20260916/mace_v1/medium/manifest.json
```

## Execute / collect

Jobs1200731/1200732 execute medium/large,20 calls each. Actual exact sbatch argv,
submission return value and manifest hash are saved in each `submission.json`.
Each uses one A5000,16 CPU,64474 MiB (one eighth of host RAM). No project
CPU/time budget; seven-day scheduler QOS request. Existing EXIT collection and
per-task receipts retain partial failures; rerunning the exact manifest reuses
only verified successful attempts. Never run two executors on one manifest.

```bash
ALQUEMIA_CURVATURE_ROOT=$(pwd -P)
sbatch --parsable --partition=gpu --nodelist=node-128-512g-8gpu-1 --job-name=alquemia_K_medium --cpus-per-task=16 --gres=gpu:1 --mem=64474M --time=7-00:00:00 --export=ALL,MACE_MIN_MEMORY_MIB=64474 --output="$ALQUEMIA_CURVATURE_ROOT/workspaces/mace_curvature_20260916/mace_v1/medium/slurm_%j.out" --error="$ALQUEMIA_CURVATURE_ROOT/workspaces/mace_curvature_20260916/mace_v1/medium/slurm_%j.err" diagnostics/mace_hybrid_20260916/run_pilot.sbatch "$ALQUEMIA_CURVATURE_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" "$ALQUEMIA_CURVATURE_ROOT/workspaces/mace_curvature_20260916/mace_v1/medium/manifest.json" native
python scripts/mace_hybrid.py collect --manifest workspaces/mace_curvature_20260916/mace_v1/medium/manifest.json --output workspaces/mace_curvature_20260916/medium_recollection.json
```

## Solvent preparation / comparison

```bash
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_curvature.py prepare-gb --collection workspaces/mace_curvature_20260916/mace_v1/medium/collection_job_1200731.json --solver-validation workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --agreement diagnostics/mace_curvature_20260916/PLAN.md --output workspaces/mace_curvature_20260916/gb_v1/medium
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_curvature.py prepare-gb --collection workspaces/mace_curvature_20260916/mace_v1/large/collection_job_1200732.json --solver-validation workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --agreement diagnostics/mace_curvature_20260916/PLAN.md --output workspaces/mace_curvature_20260916/gb_v1/large
```

The solvent jobs use the same batch wrapper/resource allocation and the literal
absolute CUDA12 solver-venv Python path. Their submitted argv/IDs and actual
collection paths are preserved in their submission records. The final report
command will name the actual completed collections once they exist.

Tests: four real-artifact preparation/algebra/failure tests passed in3.070s.
No fabricated energy/gradient fixtures or unavailable-executable replacements.

## Completed result and replay

GB jobs1200733/1200734 completed20 calls each. To regenerate the result:

```bash
python scripts/mace_curvature.py report --medium workspaces/mace_curvature_20260916/gb_v1/medium/collection_job_1200733.json --large workspaces/mace_curvature_20260916/gb_v1/large/collection_job_1200734.json --output workspaces/mace_curvature_20260916/report_replay_v1
python -m unittest discover -s tests -p test_mace_curvature.py -v
```

This reads saved artifacts only. REPORT.md/result.json/costs.json retain actual
results and receipts; no new scientific execution occurs during report replay.

Completed five-test real-fixture suite passes in3.145s, including preservation of negative curvature and null relaxation despite the passed screen.
