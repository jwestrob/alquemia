# Charge-group pilot operations

Run from the repository root. The finite scope is in [PROPOSAL.md](PROPOSAL.md).
The production baseline and all old protocols are unchanged.

## Completed execution

Jobs **1201508** (24 MACE) and **1201513** (24 OBC-II) completed. Do not submit
these manifests again. The authoritative commands are their submission.json
records under model_v3/ and solvent_v1/. Both used one A5000,16 CPUs,64474MiB.
The exact outcome and costs are in [REPORT.md](REPORT.md):1/7 directions,
30/30 numerical and3/3 grouping checks; baseline unchanged.

Two setup failures remain:1201505 failed before a forward;1201507 attempted
one forward and ran out of memory after bypassing the existing unused-grid
wrapper. V3 restores the wrappers with identical scientific inputs/settings.
All failure receipts are included in the895GPU-second total.

## Validate and collect

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_charge_groups_20260918/model_v3/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_charge_groups_20260918/model_v3/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_charge_groups.py collect --manifest workspaces/mace_charge_groups_20260918/model_v3/manifest.json
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python -m unittest discover -s tests -p test_mace_charge_groups.py -v
```

## Solvent stage — already completed; preparation reference only

The preparation command refuses missing model endpoints and does not substitute
old densities. This stage has 24 declared OBC-II tasks, with the existing solver.
The command creates a new exclusive directory and launches no computation.
The displayed solvent_v1 output already exists; do not repeat it.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_charge_groups.py prepare-gb --manifest workspaces/mace_charge_groups_20260918/model_v3/manifest.json --solver-validation workspaces/mace_gb_20260916/pilot_v2/collection_job_1200700.json --output workspaces/mace_charge_groups_20260918/solvent_v1
```

Use `software.python.path` recorded by this solvent manifest with the same
`run_pilot.sbatch` allocation; preserve its exact submission receipt. Do not use
the model Python in place of that solver environment. The executor's standard
`--task-id` selectors and existing locks/caches support partial-job recovery.

## Paired report

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_charge_groups.py report --manifest workspaces/mace_charge_groups_20260918/model_v3/manifest.json --solvent-manifest workspaces/mace_charge_groups_20260918/solvent_v1/manifest.json --native-reference workspaces/mace_global_benchmark_20260916/mace_v1/medium/collection_job_1200701.json --output workspaces/mace_charge_groups_20260918/report_review_v1
```

Omitting `--solvent-manifest` produces an explicitly incomplete vacuum-only
diagnostic; it cannot provide the declared solvent-corrected score or a pass.
Report directories are exclusive-create. Use a fresh path for any later replay.
The report retains all seven directional comparisons, numerical and
grouping failures, all attempted calculations and actual resource receipts.

## Existing preparation

`config.json` pins all seven physical/core mappings. `groups_v1/preparation.json`
records complete source groups and endpoint charges without evaluating energies.
Preparation took 11.245473 wall / 11.105136 CPU seconds. Four real tests pass
(23.170 s, zero skips), including saved native charge updates, corrupted zero
weights/charges, and preservation of existing wrapper composition. Twelve
earlier-runner tests also pass (17.125 s, zero skips). These are software/algebra
checks, not successful changed-model scientific evaluations.

`source_review_v1/preservation.json` preserves exact inspected code versions.
The pre-edit runner was recovered from its exact b3e49f5 blob and verified
against the recorded inspection hash; no scientific output was reconstructed.

## Final actual-result verification

Six real-fixture tests pass in45.957s, zero skips, including full actual report
replay, native energy/density/force identity and missing-solvent behavior.
See completed_tests.log in the workspace. additional_physical_checks_v1 contains
actual force-rotation errors; paired_comparison_v1 retains all old/new available
contrasts and the failed vacuum grouping checks. complete_cost_v1 includes all
four jobs, including failures. figure_v1 contains standalonePDF/SVG/PNG exports.
No new molecular calculation is required to replay the report or these tests.
