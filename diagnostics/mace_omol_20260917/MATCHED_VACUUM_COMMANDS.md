# Matched vacuum diagnostic commands

Protocol `native_r2scan3c_matched_vacuum_response_diagnostic_v1`; the
[plan](MATCHED_VACUUM_PLAN.md) fixes eight real centers and the comparisons.
Run from the repository root. Outputs must be fresh. No baseline/default change.

The first prepared manifest is `workspaces/mace_omol_20260917/matched_vacuum_v1/manifest.json`.
Job1200905 was submitted once; inspect it before executing any command that
could duplicate the work. Bounded archive audit searched165 actual inputs in
five named project archives:25 exact geometry matches, zero vacuum matches.
It is retained as `matched_vacuum_archive_audit_v1.json` in the parent workspace.

```bash
# Validate without executing quantum calculations.
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/matched_vacuum_v1/implementation/mace_omol_vacuum.py dry-run --manifest workspaces/mace_omol_20260917/matched_vacuum_v1/manifest.json

# Collect actual receipts; missing/invalid endpoints remain explicit.
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/matched_vacuum_v1/implementation/mace_omol_vacuum.py collect --manifest workspaces/mace_omol_20260917/matched_vacuum_v1/manifest.json --output workspaces/mace_omol_20260917/matched_vacuum_collection_replay_v1.json

# Compare actual vacuum/CPCM and BOTH previously declared learned terms.
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/mace_omol_20260917/matched_vacuum_v1/implementation/mace_omol_vacuum.py report --manifest workspaces/mace_omol_20260917/matched_vacuum_v1/manifest.json --masked workspaces/mace_omol_20260917/masked_response_report_v1/result.json --neutral workspaces/mace_omol_20260917/shared_neutral_core_report_v1/result.json --output workspaces/mace_omol_20260917/matched_vacuum_report_v1

# Real fixture tests; actual-result integration skips explicitly if absent.
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_omol_vacuum.py -v
```

Fresh preparation, if an explicit technical retry is needed:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_vacuum.py prepare --source workspaces/mace_omol_20260917/shared_neutral_core_v2/manifest.json --agreement diagnostics/mace_omol_20260917/MATCHED_VACUUM_PLAN.md --output workspaces/mace_omol_20260917/matched_vacuum_retry_v1
```

The submission wrapper takes one absolute manifest path and uses the existing
64CPU standard/memory policy, four16-rank tasks, one thread per rank. The
submitted argv and exact wrapper hash are in `matched_vacuum_v1/submission.json`.
Its trap collects results even on failure. Completed tasks can be safely resumed;
partial failures require a fresh directory and must not be overwritten. No
project budget or time limit is enforced; actual allocation costs are recorded.
