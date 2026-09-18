# Coupled donor path operations

Run from the repository root. The scientific definition is frozen in
[COUPLED_PATH_PLAN.md](COUPLED_PATH_PLAN.md). Production remains unchanged.
All compute products and actual submission argv/receipts are under
`workspaces/mace_site_response_20260918/`.

## Complete

Actual paths are `coupled_prepared_v1/preparation.json` and
`coupled_initial_v1/initial.json`. All sixteen states retain donor membership;
64 finite whole-system scalar tasks are split among five protein manifests.
Jobs 1201398/1201399 (alpha) and 1201400/1201401/1201402 (GGR 1GLG/2FVY/2FW0)
completed all64 calls. Native jobs1201407 (16DFT) and1201408 (32MACE gradients)
completed. Both frozen preflights passed. Inspect jobs before any restart.

Fresh preparation replay, with no inference:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_site_path.py prepare --projections workspaces/mace_site_response_20260918/prepared_v1/preparation.json --agreement diagnostics/mace_site_response_20260918/COUPLED_PATH_PLAN.md --output workspaces/mace_site_response_20260918/coupled_prepared_replay_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_site_path.py initial --preparation workspaces/mace_site_response_20260918/coupled_prepared_v1/preparation.json --output workspaces/mace_site_response_20260918/coupled_initial_replay_v1
```

All use the established GPU runner, one A5000 /16 CPU/64474 MiB, checkpointed
float64 analytic backend (initial tasks request only energy). Source snapshots,
software/checkpoint hashes, cache keys, failed attempts and execution locks
are retained. The GPU environment does not need the preparation-only Gemmi
package: pinned physical maps allow exact geometry replay with NumPy.

Example actual-manifest preflight and read-only collection:

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_site_response_20260918/coupled_initial_v1/ALPHA_1F6S/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_site_response_20260918/coupled_initial_v1/ALPHA_1F6S/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_site_path.py collect --manifest workspaces/mace_site_response_20260918/coupled_initial_v1/ALPHA_1F6S/manifest.json
```

## Selection and actual native checks

Only after all initial tasks complete, select the lowest of the five fixed
points, including the reused center. No change to path, radius or score rule:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_site_path.py select --initial workspaces/mace_site_response_20260918/coupled_initial_v1/initial.json --output workspaces/mace_site_response_20260918/coupled_selection_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_site_path_native.py prepare --selection workspaces/mace_site_response_20260918/coupled_selection_v1/result.json --output workspaces/mace_site_response_20260918/coupled_native_v1
```

The native preparation records explicit quantum/MACE manifests and counts.
Use `run_matched_h_quantum.sbatch ABSOLUTE_MANIFEST` and the same GPU runner;
exact submitted commands are retained beside the actual manifests. Quantum
uses64 CPU in four16-rank slots, one analytic native endpoint per selected state.
Before submitting, preflight the generated manifests with their pinned drivers.
Do not duplicate running/finished work. Existing runners retain partial failures
and reuse successful matching receipts on recovery.

Once native calculations complete:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_site_path_native.py report --preparation workspaces/mace_site_response_20260918/coupled_native_v1/preparation.json --output workspaces/mace_site_response_20260918/coupled_report_v1
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_site_path.py -v
```

The report retains all12 comparisons, all16 endpoint qualifications, both GGR
partition criteria, actual components and unavailable calibrated/free-energy
fields. The chosen point is a discrete response descriptor, not a stationary
minimum. Fresh report replays require a new output directory.

Actual selection is `coupled_selection_v1/result.json`; native manifest links
are in `coupled_native_v1/preparation.json`. The initial surrogate still gives
4/12 correct raw directions, only1GLG. Actual native results confirm4/12, with0/12 qualified comparisons.
Initial allocation cost1292GPU-s/20672core-s; native costs are additional.
Three real-fixture tests passed25.952s, none skipped in the final run.

All native jobs completed; primary result is `coupled_report_v1/result.json`.
No own job remains. See [final report](COUPLED_PATH_REPORT.md) for all failures
and actual total costs. Final four real-fixture tests pass61.767s, zero skips.
