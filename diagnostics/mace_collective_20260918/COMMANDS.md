# Collective substitution operations

Run from the repository root. Actual immutable products:
`workspaces/mace_collective_20260918/config.json`, `prepared_v1/manifest.json`,
`report_v1/result.json`, `complete_cost_v1/result.json`. Job1201453 completed all8
calls. Do not submit a duplicate executor. Source preparation failure001 is
retained; no inference failed. The existing baseline/single-site paths are unchanged.

Fresh preparation and read-only preflight/report:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_collective.py prepare --config workspaces/mace_collective_20260918/config.json --agreement diagnostics/mace_collective_20260918/PLAN.md --output workspaces/mace_collective_20260918/prepared_replay_v1
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_collective_20260918/prepared_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_collective_20260918/prepared_v1/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_collective.py report --manifest workspaces/mace_collective_20260918/prepared_v1/manifest.json --output workspaces/mace_collective_20260918/report_replay_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_collective.py -v
```

The executed command and allocation are recorded in `prepared_v1/submission.json`:
existing run_pilot.sbatch, literal pinned MACE venv Python, absolute manifest,
`native`, oneA5000/16CPU/64474MiB. Same existing execute/collect/cache/lock workflow;
no second runner. `mace_collective.py collect --manifest` reports all accepted,
missing and failed tasks with execution receipts. Explicit output directories
are exclusive-create; use new paths for report/preparation replays.
