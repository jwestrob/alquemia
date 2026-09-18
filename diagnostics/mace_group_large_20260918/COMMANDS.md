# Large grouped-model run

Job **1201533 completed** all 57 tasks. Do not submit a duplicate. The medium
PQQ job 1201524 also completed. Read PLAN.md and REPORT.md;
this is consumed development data, not a new blind benchmark.

Protocol: `intact_POLAR_large_typed31_group_vacuum_compatibility_v1`.
Manifest: `workspaces/mace_group_large_20260918/model_v1/manifest.json`.
SHA256: `abedfd09d430d9533a14325241fbbbf550943d7b5342e65c4a3910ae844c39d2`.
Exact submission arguments/allocation are in its sibling `submission.json`.

All commands run from the repository root and use fresh output paths. Existing
directories/results must remain immutable. Baseline/default are unchanged.

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_group_large_20260918/model_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_group_large_20260918/model_v1/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_group_large.py native-gate --manifest workspaces/mace_group_large_20260918/model_v1/manifest.json
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_group_large_20260918/model_v1/implementation/mace_hybrid.py collect --manifest workspaces/mace_group_large_20260918/model_v1/manifest.json --output workspaces/mace_group_large_20260918/collection_review_v1.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_group_large.py report --manifest workspaces/mace_group_large_20260918/model_v1/manifest.json --output workspaces/mace_group_large_20260918/report_review_v1
```

The first two tasks test one-group identity against actual native-large GGR
receipts. The existing executor stops before grouped tasks if that check fails.
On technical recovery it retains every attempt and reuses only accepted cache
entries for the same task/model/implementation; missing values remain null.

The report reuses the exact22fixed comparisons and source-family handling from
the medium study. `checkpoint_comparison.json` retains both margins, every gain
and regression, and the declared improvement rule (>17correct, all7original
directions retained). Numerical/native and grouping checks remain separate.
No absolute bands or response correction are provided.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_group_large.py -v
```

Prelaunch: four real source/cache/gate tests pass69.112s,zero skips; eight
medium-source/report regressions pass97.334s,zero skips. MACE-environment
preflight passed. Both actual native controls now pass: energies reproduce
exactly, maximum force difference7.084e-13eV/A and density difference3.043e-14.
Saved check: `workspaces/mace_group_large_20260918/native_gate_v1.json`.
Full report: `workspaces/mace_group_large_20260918/report_v1/`.
All 69 numerical checks pass; all three grouping checks fail. Both checkpoints
give 17/22 raw directions, with the same failures. Six complete-result/source
tests pass in 130.358 seconds, zero skipped. See VALIDATION.md and REPORT.md.
