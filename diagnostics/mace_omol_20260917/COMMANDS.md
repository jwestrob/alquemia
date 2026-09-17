# OMOL operations

Run from the repository root. All outputs are versioned and writers refuse
overwrite. Use PLAN.md for frozen scope and interpretation. Software is the
existing isolated MACE0.3.16 environment, unchanged; the new422MB checkpoint
and its official SHA256 are under workspaces/mace_omol_20260917/software_v1/.
The download receipt, release metadata and actual83-element inspection are there.
No quantum calculation, production scorer change or automatic rescore is implied.

## Prepare a new qualification replay

```bash
python scripts/mace_omol.py prepare-qualification --source-inventory workspaces/mace_canonical_20260916/audit_v2/inventory.json --software workspaces/mace_omol_20260917/software_v1/software_manifest.json --agreement diagnostics/mace_omol_20260917/PLAN.md --output workspaces/mace_omol_20260917/qualification_replay_v1
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/qualification_replay_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_20260917/qualification_replay_v1/manifest.json
```

## Execute a finite replay

The actual qualification job1200797 completed10calls and passed all9 numerical
checks. Actual benchmark job1200799 completed60newcalls and4qualified reuses. Check
its live state before any recovery; never submit concurrent executors or remove
locks. The exact submitted argv is saved beside each manifest in submission.json.
This command executes the fresh qualification replay above:

```bash
ALQUEMIA_OMOL_ROOT=$(pwd -P)
sbatch --parsable --partition=gpu --nodelist=node-128-512g-8gpu-1 --job-name=alquemia_omol_replay --cpus-per-task=16 --gres=gpu:1 --mem=64474M --time=7-00:00:00 --export=ALL,MACE_MIN_MEMORY_MIB=64474 --output="$ALQUEMIA_OMOL_ROOT/workspaces/mace_omol_20260917/qualification_replay_v1/slurm_%j.out" --error="$ALQUEMIA_OMOL_ROOT/workspaces/mace_omol_20260917/qualification_replay_v1/slurm_%j.err" diagnostics/mace_hybrid_20260916/run_pilot.sbatch "$ALQUEMIA_OMOL_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" "$ALQUEMIA_OMOL_ROOT/workspaces/mace_omol_20260917/qualification_replay_v1/manifest.json" native
```

The wrapper invokes the copied implementation and collects on exit. A valid
completed attempt is reused; partial/failed attempts remain visible. Scheduler
QOS time is not a project compute budget. All actual allocation costs count.
No host-offload or changed numerical kernel is part of this model version.

## Prepare the conditional benchmark

Only actual passing qualification can prepare the benchmark. These explicit
inputs reuse the completed qualification and original physical preparations:

```bash
python scripts/mace_omol.py prepare-benchmark --qualification workspaces/mace_omol_20260917/qualification_v1/collection_job_1200797.json --mechanics workspaces/mace_mechanics_20260916/prepared_v2/manifest.json --agreement diagnostics/mace_omol_20260917/PLAN.md --output workspaces/mace_omol_20260917/benchmark_replay_v1
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/benchmark_replay_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_20260917/benchmark_replay_v1/manifest.json
```

The fresh benchmark replay executes through the same wrapper:

```bash
ALQUEMIA_OMOL_ROOT=$(pwd -P)
sbatch --parsable --partition=gpu --nodelist=node-128-512g-8gpu-1 --job-name=alquemia_omol_bench_replay --cpus-per-task=16 --gres=gpu:1 --mem=64474M --time=7-00:00:00 --export=ALL,MACE_MIN_MEMORY_MIB=64474 --output="$ALQUEMIA_OMOL_ROOT/workspaces/mace_omol_20260917/benchmark_replay_v1/slurm_%j.out" --error="$ALQUEMIA_OMOL_ROOT/workspaces/mace_omol_20260917/benchmark_replay_v1/slurm_%j.err" diagnostics/mace_hybrid_20260916/run_pilot.sbatch "$ALQUEMIA_OMOL_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" "$ALQUEMIA_OMOL_ROOT/workspaces/mace_omol_20260917/benchmark_replay_v1/manifest.json" native
```

The actual benchmark_v1 is complete; its successful attempts will be reused
if explicitly recovered through the existing runner.

## Collect and report without inference

```bash
python scripts/mace_hybrid.py collect --manifest workspaces/mace_omol_20260917/benchmark_v1/manifest.json --output workspaces/mace_omol_20260917/benchmark_report_collection_v1.json
python scripts/mace_omol_report.py --collection workspaces/mace_omol_20260917/benchmark_report_collection_v1.json --mechanics-assessment workspaces/mace_mechanics_20260916/assessment_v1/result.json --output workspaces/mace_omol_20260917/report_replay_v1
python -m unittest discover -s tests -p test_mace_omol.py -v
```

The report verifies current receipts and preserves every denominator. The PQQ
bands use25 calibration rows only; non-PQQ directions are separate. Missing
calibration gives unavailable decisions. OMOL is a vacuum descriptor with no
predicted atomic density here: no POLAR charges, GB term, CPCM energy or old
reference is substituted. Read the completed report before any new use.
