# Coupled mechanics operations

Run from the repository root. PLAN.md and NUMERICAL_CONVENTION.md are immutable
scientific specifications. All output writers refuse overwrite. The baseline
is unchanged; these are consumed development systems with unavailable response
scores until all gates and conditional minimum validation pass.

## Preparation and audit

The executed preparation is prepared_v2. prepared_v1 stopped on an exclusive
filename collision before execution and remains preserved. Rebuild a fresh copy:

```bash
python scripts/mace_mechanics.py prepare --source-root workspaces/ggr_mechanism_20260915/stage_a_prepared_v1 --global-preparation workspaces/mace_global_benchmark_20260916/prepared_v1/preparation_manifest.json --ggr-dft workspaces/ggr_mechanism_20260915/report_v1/collection_c.json --agreement diagnostics/mace_mechanics_20260916/PLAN.md --output workspaces/mace_mechanics_20260916/preparation_replay_v1
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python scripts/mace_hybrid.py dry-run --manifest workspaces/mace_mechanics_20260916/core_v1/manifest.json
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python scripts/mace_hybrid.py dry-run --manifest workspaces/mace_mechanics_20260916/short_v1/manifest.json
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_hybrid.py dry-run --manifest workspaces/mace_mechanics_20260916/gb_v1/manifest.json
```

Preparation requires the CPU Python with gemmi. The isolated MACE Python lacks
gemmi and is deliberately used only for execution/validation paths that do not
import the molecular-preparation modules. Keep venv Python paths literal; do not
resolve their symlinks into another environment.

## Package the exact cheap inventories

These commands write fresh replay manifests; they do not launch calculations.

```bash
python scripts/mace_mechanics_run.py prepare-core --prepared workspaces/mace_mechanics_20260916/prepared_v2/manifest.json --reference-mace workspaces/mace_curvature_20260916/mace_v1/medium/collection_job_1200731.json --reference-gb workspaces/mace_curvature_20260916/gb_v1/medium/collection_job_1200733.json --short-gate workspaces/mace_short_engine_20260916/pilot_v1/collection_job_1200736.json --agreement diagnostics/mace_mechanics_20260916/PLAN.md --output workspaces/mace_mechanics_20260916/core_replay_v1
python scripts/mace_mechanics_run.py prepare-short --prepared workspaces/mace_mechanics_20260916/prepared_v2/manifest.json --short-gate workspaces/mace_short_engine_20260916/pilot_v1/collection_job_1200736.json --agreement diagnostics/mace_mechanics_20260916/PLAN.md --output workspaces/mace_mechanics_20260916/short_replay_v1
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_mechanics_run.py prepare-gb --collection workspaces/mace_mechanics_20260916/core_v1/collection_job_1200749.json --solver-validation workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --agreement diagnostics/mace_mechanics_20260916/PLAN.md --output workspaces/mace_mechanics_20260916/gb_replay_v1
```

## Execution and recovery

Actual sbatch argv, IDs and input hashes are in each submission.json. Quantum
job1200743 runs36 new DFT tasks through affordable_workflow, four16-rank workers
on64CPU. MACE core1200749 runs116 new calls; short1200750 runs106 calls; solvent
1200767 runs116 calls. Each GPU job uses one A5000,16CPU and64474MiB. Cached20GGR
core/solvent calls are independently pinned to their existing receipts. No new
DFT is hidden inside a cheap task. No project CPU/time stopping budget exists.

Existing wrappers are diagnostics/mace_hybrid_20260916/run_pilot.sbatch and
run_quantum.sbatch in this directory. Never resubmit a live manifest. The MACE
wrapper uses the manifest's implementation snapshot and writes a collection on
exit. Valid completed attempts are reused; failures stay visible. Do not remove
locks or mutate a live snapshot. Quantum partial attempts require the existing
runner's explicit fresh retry mechanism, not overwritten outputs.

Read-only collection and real-fixture tests:

```bash
python scripts/mace_hybrid.py collect --manifest workspaces/mace_mechanics_20260916/core_v1/manifest.json --output workspaces/mace_mechanics_20260916/core_recollection_v1.json
python scripts/mace_hybrid.py collect --manifest workspaces/mace_mechanics_20260916/short_v1/manifest.json --output workspaces/mace_mechanics_20260916/short_recollection_v1.json
python scripts/mace_hybrid.py collect --manifest workspaces/mace_mechanics_20260916/gb_v1/manifest.json --output workspaces/mace_mechanics_20260916/gb_recollection_v1.json
python scripts/mace_mechanics.py collect-dft --manifest workspaces/mace_mechanics_20260916/prepared_v2/quantum/manifest.json --output workspaces/mace_mechanics_20260916/quantum_recollection_v1.json
python -m unittest discover -s tests -p 'test_mace_mechanics*.py' -v
```

## Compare and report

After all four exact collections exist and pass execution validation:

```bash
python scripts/mace_mechanics_assess.py --prepared workspaces/mace_mechanics_20260916/prepared_v2/manifest.json --core workspaces/mace_mechanics_20260916/core_v1/collection_job_1200749.json --gb workspaces/mace_mechanics_20260916/gb_v1/collection_job_1200767.json --short workspaces/mace_mechanics_20260916/short_v1/collection_job_1200750.json --dft workspaces/mace_mechanics_20260916/prepared_v2/quantum/collection_job_1200743.json --normal-manifest workspaces/ggr_mechanism_20260915/stage_a_tasks_v1/manifest.json --convention diagnostics/mace_mechanics_20260916/NUMERICAL_CONVENTION.md --output workspaces/mace_mechanics_20260916/assessment_replay_v1
```

The assessment checks paired states, energy/gradient mappings, DFT local
curvature, grid refinement, normal/tight SCF bridges, eigenvalues and trust
limits. It reports predicted minima and all failures. It does not turn a
predicted minimum into a validated relaxation score. Only eligible fixed minima
can advance to the already planned independent DFT validation (up to8 tasks).


## Fixed-minimum validation (three actual eligible Ca endpoints)

The initial assessment is assessment_v1/result.json. All prerequisites pass,
but only alpha1F6S Ca and both GGR Ca representations are eligible. Preparation
selects every eligible endpoint and rejects status-only promotion of an
unsupported one. Actual jobs1200771 (3DFT) and1200772 (6short calls) completed.
Their submission records contain the exact commands and prediction hashes.

```bash
python scripts/mace_mechanics_minimum.py prepare --assessment workspaces/mace_mechanics_20260916/assessment_v1/result.json --output workspaces/mace_mechanics_20260916/minimum_replay_v1
python scripts/mace_mechanics_minimum.py report --prepared workspaces/mace_mechanics_20260916/minimum_v1/preparation.json --dft workspaces/mace_mechanics_20260916/minimum_v1/quantum/collection_job_1200771.json --short workspaces/mace_mechanics_20260916/minimum_v1/short/collection_job_1200772.json --output workspaces/mace_mechanics_20260916/minimum_report_replay_v1
python -m unittest discover -s tests -p test_mace_mechanics_minimum.py -v
```

No DFT gradient was calculated at the predicted minimum: the validation concerns
its actual energy change, not proof that the point is an exact DFT stationary
minimum. All three energy checks pass; missing La corrections remain null.
