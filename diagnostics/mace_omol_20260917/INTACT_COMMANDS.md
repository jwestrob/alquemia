# Intact-chain OMOL operations

Read INTACT_CHAIN_PLAN.md for the fixed inputs and interpretation. Production
is unchanged. Energy-only output explicitly has no forces. All task manifests
use the existing immutable runner, receipts and typed cache validation.

## Core equivalence bridge

The four-task manifest was prepared and passed the copied runner dry-run.
Job1200807 uses `intact_core_v1`; its exact sbatch argv is saved in
`intact_core_v1/submission.json`. These commands inspect without new inference:

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/intact_core_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_20260917/intact_core_v1/manifest.json
python -m unittest discover -s tests -p test_mace_omol_intact.py -v
```

## Intact numerical qualification

After the actual core collection passes, prepare the14-task ALPHA_1F6S gate:

```bash
python scripts/mace_omol_intact.py --collection workspaces/mace_omol_20260917/benchmark_v1/collection_job_1200799.json --preparation workspaces/mace_global_benchmark_20260916/prepared_v1/preparation_manifest.json --agreement diagnostics/mace_omol_20260917/INTACT_CHAIN_PLAN.md --core-qualification workspaces/mace_omol_20260917/intact_core_v1/collection_job_1200807.json --output workspaces/mace_omol_20260917/intact_qualification_v1
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/intact_qualification_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_20260917/intact_qualification_v1/manifest.json
```

Preparation refuses overwrite. Do not repeat preparation for an existing
manifest; inspect its receipt and live job before recovery. Core bridge success
is an engineering check, not predictive validation. The next16-task comparison
is conditional on the14-task numerical gate, irrespective of score direction.

Qualification manifest SHA256:
`b1609f0b81c73592881102ca123bce39f78d529ba23b64373d26d9e83d491488`.
Preparation ran the same validator as dry-run and passed. Job1200808 executes
this manifest; exact allocation/submission arguments are in its `submission.json`.
The executor repeats validation before inference and collects on exit.

```bash
squeue -j 1200808 -o '%.12i %.24j %.10T %.10M %.24R'
python scripts/mace_hybrid.py collect --manifest workspaces/mace_omol_20260917/intact_qualification_v1/manifest.json --output workspaces/mace_omol_20260917/intact_qualification_manual_collection_v1.json
```

Manual collection computes no new scientific outputs. It preserves missing
endpoints as unavailable; it cannot make an incomplete numerical gate pass.


## Memory recovery

The first native full call in1200808 exhausted A5000 GPU memory; no energy was
produced. Job1200809 is the unchanged-manifest H200 recovery (28 CPUs,
200000 MiB host entitlement). Its exact argv is in
`intact_qualification_v1/submission_h200_recovery.json`; do not duplicate it.
The failed attempt remains part of the manifest's execution history.

Read EXACT_EDGE_PLAN.md for the independently declared memory qualification.
The adapter retains every edge and original nonlinear message aggregation.
Prepare its eight-call core test with these actual input artifacts:

```bash
python scripts/mace_omol_edge_run.py --native-core-collection workspaces/mace_omol_20260917/intact_core_v1/collection_job_1200807.json --intact-manifest workspaces/mace_omol_20260917/intact_qualification_v1/manifest.json --agreement diagnostics/mace_omol_20260917/EXACT_EDGE_PLAN.md --output workspaces/mace_omol_20260917/edge_core_replay_v1
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/edge_core_replay_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_20260917/edge_core_replay_v1/manifest.json
python -m unittest discover -s tests -p test_mace_omol_edges.py -v
```

A fresh manifest uses the existing run_pilot.sbatch executor. Actual submitted
jobs and their exact arguments will be recorded beside their own manifests.
A native cache cannot satisfy an adapted task, or conversely. Unavailable native
intact equivalence is not a pass. The adapter currently supports energy only.


## Completed batching checks and native CPU reference

Batched core1200810 passed20 checks; intact A5000 qualification1200811 passed31.
Native CPU core1200812 passed10 checks. The14-task native CPU intact reference
is job1200814; exact allocation args are in `cpu_intact_v1/submission.json`.
It uses the same runner through `run_cpu.sbatch`,64 CPUs and128GiB host RAM.
Read NATIVE_CPU_PLAN.md and INTACT_ENGINEERING_REPORT.md before interpretation.

To prepare an explicit fresh CPU reference replay after its core gate:

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python scripts/mace_omol_edge_run.py --backend cpu --native-core-collection workspaces/mace_omol_20260917/intact_core_v1/collection_job_1200807.json --intact-manifest workspaces/mace_omol_20260917/intact_qualification_v1/manifest.json --agreement diagnostics/mace_omol_20260917/NATIVE_CPU_PLAN.md --backend-qualification workspaces/mace_omol_20260917/cpu_core_v1/collection_job_1200812.json --output workspaces/mace_omol_20260917/cpu_intact_replay_v1
```

After actual CPU completion, verify full native/batched equivalence without
new inference:

```bash
python scripts/mace_omol_edge_report.py --native-collection workspaces/mace_omol_20260917/cpu_intact_v1/collection_job_1200814.json --edge-collection workspaces/mace_omol_20260917/edge_intact_v1/collection_job_1200811.json --output workspaces/mace_omol_20260917/edge_equivalence_v1
```

Only a passing actual equivalence report permits this16-task continuation:

```bash
python scripts/mace_omol_intact.py --collection workspaces/mace_omol_20260917/benchmark_v1/collection_job_1200799.json --preparation workspaces/mace_global_benchmark_20260916/prepared_v1/preparation_manifest.json --agreement diagnostics/mace_omol_20260917/INTACT_BATCHED_CONTINUATION.md --core-qualification workspaces/mace_omol_20260917/intact_core_v1/collection_job_1200807.json --full-qualification workspaces/mace_omol_20260917/cpu_intact_v1/collection_job_1200814.json --edge-equivalence workspaces/mace_omol_20260917/edge_equivalence_v1/result.json --output workspaces/mace_omol_20260917/intact_benchmark_v1
```

The continuation reuses four verified batched ALPHA_1F6S endpoints. Its other
16tasks retain every frozen physical state and the original three relative
criteria. No predictive result or calibrated class is available yet.
