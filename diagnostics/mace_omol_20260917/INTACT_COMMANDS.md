# Intact-chain OMOL operations

Read INTACT_CHAIN_PLAN.md for the fixed inputs and interpretation. Production
is unchanged. Energy-only output explicitly has no forces. All task manifests
use the existing immutable runner, receipts and typed cache validation.

## Current checkpoint — 2026-09-17

The five-chain benchmark completed in jobs1200815/1200816; all three relative
criteria pass. Read INTACT_REPORT.md and INTACT_RESULT.json. The original
intact_report_v1 is immutable. A fresh report replay (no new inference) is:

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/intact_reporting_source_v1/implementation/mace_omol_intact_report.py --collection workspaces/mace_omol_20260917/intact_benchmark_v1/collection_job_1200816.json --output workspaces/mace_omol_20260917/intact_report_replay_v1
```

Product batching core1200817 passed20checks exactly. Full qualification1200818
completed14 calls and31 checks; all24 native CPU/product GPU comparisons agree
exactly. Its manifestSHA is
`6ece3a194ab1e88cbefab44d62c3fee39cde7bdf98e6594caf4fcc5952f4f2eb`.
The completed comparison used this frozen source (output already exists):

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/product_reporting_source_v1/implementation/mace_omol_edge_report.py --native-collection workspaces/mace_omol_20260917/cpu_intact_v1/collection_job_1200814.json --edge-collection workspaces/mace_omol_20260917/product_intact_v1/collection_job_1200818.json --output workspaces/mace_omol_20260917/product_equivalence_v1
```

The canonical experiment declared in INTACT_CANONICAL_PLAN.md is running as
job1200819:104 new endpoints and8 exact earlier crystal reuses. Its immutable
manifest is `intact_panel_run_v1/manifest.json`, SHA256
`642d9717d69a8ee236e1fbdc280fa9fbed2d169e2c073c29efac7f63a72cca62`.
Four calibration charges exceed the reported training range.

**Mandatory preparation audit:** 1KB0 contains two false peptide connections
across missing structure. Its four raw endpoints are invalid preparation
diagnostics, even if inference succeeds. Read
INTACT_CANONICAL_INTEGRITY_ADDENDUM.md. Retain the unavailable case in the
three-transfer denominator. Use the new audited reporter, never the original
reporter frozen with the running manifest:

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/intact_panel_reporting_source_v2/implementation/mace_omol_panel_report.py --collection workspaces/mace_omol_20260917/intact_panel_run_v1/collection_job_1200819.json --preparation-audit workspaces/mace_omol_20260917/intact_panel_integrity_v1/result.json --output workspaces/mace_omol_20260917/intact_panel_report_v1
```

Run after the collection exists; this computes no new endpoints. The v2
preparation under `intact_panel_prepared_v2/` rejects 1KB0 before template
matching and retains27 supported cases. It does not replace the running job's
immutable v1 inputs. Root preparation commands also reject the invalid legacy
preparation when building a new run. Eight real-fixture tests passed.

Earlier readiness/template success was insufficient to detect these gaps.
Do not rerun the historical preparation/submission commands below into existing
outputs, or launch duplicate jobs. They document how the earlier stages ran.

## Historical core equivalence bridge

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


Full native equivalence passed all24 comparisons exactly. Native CPU1200814 is
complete; benchmark1200815 is running16newcalls/four qualified batched reuses.
Manifest SHA256:
`2451a00c17b27b0ba9e850312fd7684d38fd255f64e605b3402422a22864efb4`.
Exact submission argv is in `intact_benchmark_v1/submission.json`.
After collection, use the frozen reporting source (do not edit it while running):

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/intact_reporting_source_v1/implementation/mace_omol_intact_report.py --collection workspaces/mace_omol_20260917/intact_benchmark_v1/collection_job_1200815.json --output workspaces/mace_omol_20260917/intact_report_v1
python -m unittest discover -s tests -p test_mace_file_checks.py -v
python -m unittest discover -s tests -p test_mace_omol_edges.py -v
```

The earlier live report-source pin limitation and independent exact replay are
recorded in INTACT_REPORT_SOURCE_NOTE.md. No inference result changed.
