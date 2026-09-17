# Whole-protein development runbook

Use repository-relative paths below from the repository root. Existing outputs
are immutable; collection/preparation refuse to overwrite. Baseline defaults
are unchanged. Full scope and exclusions are in PLAN.md.

## Prepared and running artifacts

Physical inputs: `workspaces/mace_global_benchmark_20260916/prepared_v1/`.
Medium manifest has14 tasks; large has10. Both passed dry-run. Jobs1200701 and
1200702 execute these finite manifests with one A5000,16CPU,64474MiB each.

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python scripts/mace_hybrid.py dry-run --manifest workspaces/mace_global_benchmark_20260916/mace_v1/medium/manifest.json
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python scripts/mace_hybrid.py dry-run --manifest workspaces/mace_global_benchmark_20260916/mace_v1/large/manifest.json
```

The batch EXIT handler writes `collection_job_JOBID.json`, retaining partial
failures. The existing runner resumes successful tasks from verified receipts.
Do not start another executor while its lock is held.

## Prepare and run the declared solvent stage after MACE completion

```bash
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_global_benchmark.py prepare-gb --collection workspaces/mace_global_benchmark_20260916/mace_v1/medium/collection_job_1200701.json --solver-validation workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --agreement diagnostics/mace_global_benchmark_20260916/PLAN.md --output workspaces/mace_global_benchmark_20260916/gb_v1/medium
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_global_benchmark.py prepare-gb --collection workspaces/mace_global_benchmark_20260916/mace_v1/large/collection_job_1200702.json --solver-validation workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --agreement diagnostics/mace_global_benchmark_20260916/PLAN.md --output workspaces/mace_global_benchmark_20260916/gb_v1/large
```

Use the existing `diagnostics/mace_hybrid_20260916/run_pilot.sbatch` with explicit
Python and absolute manifest paths, memory mode `native`. The CUDA12 Python is
`workspaces/mace_gb_20260916/software_v2/venv/bin/python`; retain the same A5000
allocation flags as above and `MACE_MIN_MEMORY_MIB=64474`. No MACE/DFT rerun is
part of this solvent stage. No score is filled with zero on failure.

## Reproduce preparation in a fresh directory

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_global_prepare.py --audit diagnostics/global_electrostatic_20260916/ACCURACY_INPUTS.json --plan diagnostics/mace_global_benchmark_20260916/PLAN.md --output workspaces/mace_global_benchmark_20260916/preparation_recheck
```

This is deterministic geometry/charge/topology preparation, not optimization.
It keeps every source heavy atom, completes only approved missing terminalOXT,
fixes declared protein/water H lengths and records every excluded source residue.
PQQ microstate and coordinates remain frozen. Unsupported inputs remain explicit.

## Completed collections and paired report

Both solvent jobs completed: medium1200711, large1200712. The direct model
failed all three expected-order checks for both checkpoints; see REPORT.md.
No production score was changed.

```bash
workspaces/mace_gb_20260916/software_v2/venv/bin/python scripts/mace_global_benchmark.py compare --medium-gb workspaces/mace_global_benchmark_20260916/gb_v1/medium/collection_job_1200711.json --large-gb workspaces/mace_global_benchmark_20260916/gb_v1/large/collection_job_1200712.json --output workspaces/mace_global_benchmark_20260916/comparison_recheck
```

The stored XYZ comment inherited from the original1H4I pilot incorrectly names
1H4I for every case; the actual identity/geometry is correctly pinned in each
manifest. Those immutable bytes are preserved. The writer now uses a generic
comment for future artifacts; coordinates are unaffected. The completed exact
preparation replay used the original writer and matched all XYZ hashes.
