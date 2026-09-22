# Remaining 26 angular proposals: prepared operations

The 52-task manifest is ready at
`workspaces/adaptive_accommodation_20260922/remaining26_proposals_v1/manifest.json`.
No optimization or molecular call was run by this preparation. Root owns the
execution decision after the separate numerical check. Settings, model and
physical constraints are unchanged from the completed four-case pilot.

## Validate the exact prepared snapshot

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export ADAPTIVE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$ADAPTIVE_PY" \
  workspaces/adaptive_accommodation_20260922/remaining26_proposals_v1/implementation/adaptive_angular_proposals.py \
  dry-run \
  --manifest workspaces/adaptive_accommodation_20260922/remaining26_proposals_v1/manifest.json
```

Expected status: `prepared_dry_run_pass`, 26 cases, 52 tasks, zero molecular calls.
The actual prelaunch `collection_unrun_v1.json` retains all 52 candidates as
unavailable and every new solvent/composite/score field null. It is not a failed
scientific run and must not be overwritten after a later execution.

## Reproduce preparation in a new directory

```bash
"$ADAPTIVE_PY" scripts/adaptive_angular_proposals.py prepare \
  --source workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json \
  --diagnostic workspaces/adaptive_accommodation_20260922/force_projection_v1.json \
  --agreement diagnostics/adaptive_accommodation_20260922/CONTINUATION_PREPARATION_PLAN.md \
  --population remaining26 \
  --output workspaces/adaptive_accommodation_20260922/remaining26_replay_v1
```

Omitting `--population` retains the original `pilot4` behavior. Explicit
`--population pilot4` is equivalent. This option selects a fixed source set,
not different chemistry or a different optimizer. Both populations retain the
same runtime prerequisite: the completed **original four-case** common pool at
`workspaces/nikasha_shared_pool_20260922/pilot_v2/after_solvent_0_1209845.json`.
The pending/adaptive-expanded pool is not substituted for that prerequisite.

## Existing allocation wrapper, held for root's execution decision

The existing wrapper accepts this manifest without modification. The following
is an explicit future submission command, **not executed by this task**:

```bash
ADAPTIVE_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
sbatch --parsable --chdir="$ADAPTIVE_ROOT" \
  --output="$ADAPTIVE_ROOT/diagnostics/adaptive_accommodation_20260922/remaining26_%j.out" \
  --error="$ADAPTIVE_ROOT/diagnostics/adaptive_accommodation_20260922/remaining26_%j.err" \
  "$ADAPTIVE_ROOT/diagnostics/adaptive_accommodation_20260922/run_proposals.sbatch" \
  "$ADAPTIVE_ROOT/workspaces/adaptive_accommodation_20260922/remaining26_proposals_v1/manifest.json" \
  "$ADAPTIVE_ROOT/workspaces/nikasha_shared_pool_20260922/pilot_v2/after_solvent_0_1209845.json"
```

It uses the existing one-H200 / 32-CPU / 200000-MiB host-memory recipe and immutable
implementation snapshot, then collects actual endpoint statuses. These are 52
optimizer starts, not a claim of 52 energy/force calls. No GFN2, DFT, calibration,
common-pool selection or production update occurs inside this adapter.

## Fixture tests

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  -m unittest discover -s tests -p 'test_adaptive_angular_*.py' -v
```

The new continuation tests use all real remaining maps and retained original
results; finite differences check only analytic geometry constraints. No energy
engine or optimizer success is mocked.
