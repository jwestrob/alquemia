# Runnable operations

**Current continuation:** Jacob removed compute budgets and time limits; see
[CONTINUATION_AGREEMENT.md](CONTINUATION_AGREEMENT.md). The original pilot and
two serial solver attempts are complete historical executions. Job 1198968
executes the remaining 102 independent charging solves concurrently and reuses
the completed state. These are the same frozen 18 scientific states.

```bash
squeue -j 1198968
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/affordable_solver.py dry-run-completion \
  --manifest workspaces/affordable_challenger_20260915/solver_completion/completion_manifest.json
```

The active script is [run_solver_completion.sbatch](run_solver_completion.sbatch).
It invokes `execute-completion --manifest` inside Slurm. Repeating a successful
charging task requires matching executable, input and output receipts before
cache reuse; partial outputs are retained. Every remaining task executes
regardless of accumulated cost. Concurrency is limited only by task count,
allocated CPUs and available memory. Do not resubmit while the job is running.

Upon completion, read `workspaces/affordable_challenger_20260915/solver_completion/REPORT.md`.
To regenerate a comparison, pass the explicit result path:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/affordable_compare.py --root . \
  --solver-result workspaces/affordable_challenger_20260915/solver_completion/solver_result.json \
  --output diagnostics/affordable_challenger_20260915/comparison_completed_checks
```

## Historical initial-pilot operations

The sections below describe the original manifests. Their budgets and endpoint
caps are historical metadata and are no longer enforced. The original solver
launch scripts require their preserved implementation snapshots; the current
solver CLI uses `prepare-completion`, `dry-run-completion`, `execute-completion`.

Run on biotite. All output operations create new files/directories and refuse
overwrites. The two pilot jobs are **already submitted**; monitor them rather
than submitting duplicate jobs. Existing production commands are unchanged.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ALQUEMIA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Audit, dry-run, collect, compare and report

```bash
"$ALQUEMIA_PY" scripts/affordable_audit.py --root . \
  --output diagnostics/affordable_challenger_20260915/live_audit_refresh.json
"$ALQUEMIA_PY" scripts/affordable_workflow.py dry-run \
  --manifest workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json
squeue -j 1198934,1198939
"$ALQUEMIA_PY" scripts/affordable_workflow.py collect \
  --manifest workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json \
  --output diagnostics/affordable_challenger_20260915/collection_review.json
"$ALQUEMIA_PY" scripts/affordable_compare.py --root . \
  --output diagnostics/affordable_challenger_20260915/comparison_after_jobs
```

`compare` writes both JSON records and `REPORT.md`; no fitting or new compute.
Collecting before completion retains unavailable statuses. Job execution
automatically preserves source copies and records actual hashes, receipts and
allocation usage. The submitted solver also writes its final comparison/report.

## Reproduce preparation without submitting calculations

```bash
"$ALQUEMIA_PY" scripts/affordable_peptide.py \
  --baseline-manifest diagnostics/nonpqq_direct_site_benchmark_20260915/prepared/ggr_1glg/GGR/ggr_1glg_GGR_carve_manifest.json \
  --topology /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/lib/python3.11/site-packages/openmm/app/data/amber19/protein.ff19SB.xml \
  --output workspaces/affordable_challenger_20260915/ggr_command_reproduction
"$ALQUEMIA_PY" scripts/affordable_workflow.py prepare --root . \
  --output workspaces/affordable_challenger_20260915/pilot_command_reproduction
"$ALQUEMIA_PY" scripts/affordable_state.py \
  --pilot-manifest workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json \
  --output workspaces/affordable_challenger_20260915/environment_command_reproduction
```

The remaining five exact repair sources are in `repair_inventory_verified.json`.
The current verified workspace is `verified_repairs/`. Environment preparation
records unsupported states and measured process cost; it does not synthesize
missing charges or run dynamics. Do not execute a reproduced pilot as extra
development endpoints: the eight-endpoint cap is campaign-wide.

## Bounded execution and recovery

The actual submitted launch scripts are [run_pilot.sbatch](run_pilot.sbatch)
and [run_solver.sbatch](run_solver.sbatch); submission receipts live in the
task workspace. They specify existing CPU partitions and MPI policy, invoke
the existing manifested runner, and enforce admission budgets. They must run
inside Slurm, not on the login node. No cluster or production defaults changed.

The underlying bounded operation, already scheduled by the first script, is:

```bash
"$ALQUEMIA_PY" scripts/affordable_workflow.py execute \
  --manifest workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json
```

Successful entries require matching inputs, implementation hashes and execution
receipts for cache reuse. Partial attempts are preserved. `prepare-retry
--manifest PATH --task ID` accepts only an actually failed original endpoint
and creates a fresh identical attempt sharing the original budget ledger;
no retry is currently needed or prepared. An interrupted execution with unknown
allocated cost blocks more admissions until accounting is reconciled. Retry
results require explicit lineage-aware collection; there is no silent choice
of a favorable attempt. Low-level failed runs likewise require explicit review
and shared-budget reconciliation, not a fresh unaccounted solver directory.

## Tests and response interface

```bash
"$ALQUEMIA_PY" -m unittest discover -s tests -p test_affordable_development.py -v
"$ALQUEMIA_PY" scripts/affordable_response.py --help
```

The response CLI accepts explicit `--engrad`, `--orca-output`, `--input`,
`--xyz`, optional `--repair-manifest`, and `--output`. The paired diagnostic
accepts `--la-gradient`, `--ca-gradient`, `--output`, requiring matching physical
source mappings. No matching real analytic-gradient artifacts exist for this
pilot, so no numerical extraction command is represented as executed.
Mechanical and entropy values remain null: `response_model_not_validated`.

## Provenance / dependencies

Use `IMPLEMENTATION_REVIEW.json` for exact live source snapshot and installed
Python/package versions. The ORCA binary and runner hashes are pinned in the
pilot manifest. The APBS official archive is
`workspaces/affordable_challenger_20260915/software/APBS-3.4.1.Linux.zip`,
SHA-256 `750f6a2df7b5a82b69be5c4cb192115c02fa261d0ff52cf1d9632ed1eda8b4f0`.
Its isolated executable is under `software/apbs-3.4.1/APBS-3.4.1.Linux/bin/apbs`.
No shared environment installation was modified.
