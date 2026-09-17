# Exact short engine operations

Run from the repository root. All outputs refuse overwrite. Read PLAN.md and
REPORT.md first; the executed34-task experiment is complete, not queued work.

## Prepare a fresh replay / audit

```bash
python scripts/mace_short_engine.py --full workspaces/mace_global_benchmark_20260916/mace_v1/medium/collection_job_1200701.json --core workspaces/mace_curvature_20260916/mace_v1/medium/collection_job_1200731.json --agreement diagnostics/mace_short_engine_20260916/PLAN.md --output workspaces/mace_short_engine_20260916/replay_v1
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python scripts/mace_hybrid.py dry-run --manifest workspaces/mace_short_engine_20260916/pilot_v1/manifest.json
```

## Execute / recover / collect

The existing MACE runner executes the finite manifest and caches only verified
successful attempts. The submitted argv is in pilot_v1/submission.json; resource
allocation was one A5000,16CPU,64474MiB. Never run concurrent executors on one
manifest or remove another executor's lock. For the fresh replay prepared above:

```bash
ALQUEMIA_SHORT_ROOT=$(pwd -P)
sbatch --parsable --partition=gpu --nodelist=node-128-512g-8gpu-1 --job-name=alquemia_short_replay --cpus-per-task=16 --gres=gpu:1 --mem=64474M --time=7-00:00:00 --export=ALL,MACE_MIN_MEMORY_MIB=64474 --output="$ALQUEMIA_SHORT_ROOT/workspaces/mace_short_engine_20260916/replay_v1/slurm_%j.out" --error="$ALQUEMIA_SHORT_ROOT/workspaces/mace_short_engine_20260916/replay_v1/slurm_%j.err" diagnostics/mace_hybrid_20260916/run_pilot.sbatch "$ALQUEMIA_SHORT_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" "$ALQUEMIA_SHORT_ROOT/workspaces/mace_short_engine_20260916/replay_v1/manifest.json" native
python scripts/mace_hybrid.py collect --manifest workspaces/mace_short_engine_20260916/pilot_v1/manifest.json --output workspaces/mace_short_engine_20260916/recollection_v1.json
python -m unittest discover -s tests -p test_mace_short_engine.py -v
```

Collect runs the component, rotation and physical-direction comparisons on actual
saved results; it launches no inference. The historical collection is
pilot_v1/collection_job_1200736.json. The seven-day request follows scheduler QOS;
there is no project CPU-time stopping budget. Preparation/task parameters are
recorded in the manifest, not inferred from shell state.
