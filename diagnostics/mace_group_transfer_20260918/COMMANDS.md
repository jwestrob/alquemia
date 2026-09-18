# Typed-contact group transfer operations

Run from the repository root. The immutable scope is in PLAN.md; sources in
INPUTS.json. This is opt-in research; production scoring/reference unchanged.

## Completed execution

**Job1201517 completed all55evaluations. Do not resubmit it.** Exact command:
workspaces/mace_group_transfer_20260918/model_v3/submission.json. Finalreport
in report_v1/result.json:17/22 rawdirections,61/61numericchecks,0/3groupingchecks.
[Interpretation and actualcosts](REPORT.md). No ownlivejobs.

## Current preparation

All34cases are prepared. The finite55-task manifest is model_v3/manifest.json
under workspaces/mace_group_transfer_20260918. Scientific settings and all55XYZ
hashes match V1/V2; preflights fixed a missing gemmi import by loading exact
pure donor definitions and a4e-16A norm-rounding mismatch. No model calculation
ran in those failed preflights. Source_review_v1 and failure receipts preserve
both. Prepared_v2 snapshots its source instead of relying on mutable code paths.

## Validate / collect

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_group_transfer_20260918/model_v3/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_group_transfer_20260918/model_v3/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_group_transfer.py collect --manifest workspaces/mace_group_transfer_20260918/model_v3/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_group_transfer.py -v
```

The first two commands read existing state; tests use only real pinned sources
and explicitly corrupted copies for failure handling. Shared all-Ca references
are admitted only after physical atom/group equality across all sites.

## Execute

Use the existing diagnostics/mace_hybrid_20260916/run_pilot.sbatch with the
pinned model Python and exact manifest above. One A5000,16CPUs,64474MiB. The
submission receipt records the executable command and actual job ID; inspect
that job before any retry. Existing locks, accepted-result caching and optional
--task-id selectors support recovery without replacing failed attempts.
A successful baseline or older group-model cache cannot satisfy these tasks.

## Report

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_group_transfer.py report --manifest workspaces/mace_group_transfer_20260918/model_v3/manifest.json --output workspaces/mace_group_transfer_20260918/report_review_v1
```

This writes an exclusive result directory. All22directional comparisons,three
GGR grouping controls,all55charge checks,three rigid/permutation controls and
ordered aequorin scores remain visible. Missing endpoints remain unavailable;
no solvent/baseline fallback or fitted threshold. Repeat with a fresh output
path. No scientific calculation is triggered by report generation.

## Input/model preparation already performed

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_group_transfer.py prepare-inputs --config workspaces/mace_group_transfer_20260918/config.json --output workspaces/mace_group_transfer_20260918/prepared_v2
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_group_transfer.py prepare --preparation workspaces/mace_group_transfer_20260918/prepared_v2/preparation.json --output workspaces/mace_group_transfer_20260918/model_v3
```

Those output directories exist: do not rerun into them. Their exact config and
snapshots support a fresh replay with another explicit output directory.

Six real source tests pass60.417s,zero skips. Six previous-group regression tests
pass45.119s,zero skips. The MACE-environment frozen preflight passes. Primary
scientific comparison and actual-result integration tests are now complete.

## Completed actual-result tests

Eighttests pass96.968s,zero skips. Includes actualfullreport replay,partial-output
unavailability,means,orderedsitevectors,sign/unit algebra andreal source-mapping
checks. See completed_tests.log; partial_report_v1/result.json preserves a real
24/55endpoint snapshot with missing scores unavailable. No fictitious success.
Completecost_v1 includesactualGPU/CPUreceipts andboth preflight/preparation
attempts. Comparison_v1 providesold/new supportingmargins andcomponentaudit;
figure_v1 containsPDF/SVG/PNG. Report replay addsno scientific evaluation.
