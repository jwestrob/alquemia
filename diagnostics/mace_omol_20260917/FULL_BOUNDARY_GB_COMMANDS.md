# Full physical boundary solvent challenger: operations

From the repository root. All new scientific outputs must use fresh paths.
Job1200975 was submitted once; inspect its live state and receipts before any
execution. This is an uncalibrated research descriptor, with no production change.
[Plan](FULL_BOUNDARY_GB_PLAN.md) fixes physical inputs, charge boundary and48tasks.

Validate/collect without additional solver execution:

```bash
workspaces/mace_gb_20260916/software_v2/venv/bin/python workspaces/mace_omol_20260917/full_boundary_GB_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_20260917/full_boundary_GB_v1/manifest.json

workspaces/mace_gb_20260916/software_v2/venv/bin/python workspaces/mace_omol_20260917/full_boundary_GB_v1/implementation/mace_omol_solvent.py collect --manifest workspaces/mace_omol_20260917/full_boundary_GB_v1/manifest.json --output workspaces/mace_omol_20260917/full_boundary_GB_collection_replay_v1.json

/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_omol_solvent.py -v
```

Fresh preparation without cluster submission:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_solvent.py prepare --charges workspaces/mace_omol_20260917/normalized_charge_report_v1/result.json --hybrid workspaces/mace_omol_20260917/matched_H_report_v1/result.json --software workspaces/mace_gb_20260916/pilot_v2/software.json --agreement diagnostics/mace_omol_20260917/FULL_BOUNDARY_GB_PLAN.md --output workspaces/mace_omol_20260917/full_boundary_GB_preparation_replay_v1
```

The original submission argv are retained in
`workspaces/mace_omol_20260917/full_boundary_GB_v1/submission.json`.
Use the existing `run_pilot.sbatch`, pinned GB CUDA12 Python, A5000 allocation,
16CPU/64474MiB and recorded memory override. The shared runner executes finite
tasks, reuses verified successes, keeps failures and collects on exit. It never
substitutes a baseline result for a failed solvent task. No application compute
budget is imposed. Never rerun completed solver tasks merely to regenerate a
report. Inspect any failed receipt before an explicit technical retry.

Source charge fitting and cap projection are already complete. New state files
retain the original ff19SB charges, excluded QM support, local recipients and
charge increments, exact full coordinates, endpoint charge vectors, formal
ledgers and evidence strata. Derived component charges are mathematical sources
on the full cavity, not different protonation/electronic states. Solvent forces
hold charges fixed and are unavailable as combined-model gradients.
