# Approved parallel pilots: operations

All commands run from the repository root, with the pinned CPU interpreter:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
NIKASHA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
```

Each branch saves its specification, finite manifest, executor snapshot and
submission receipts under its own workspace. The prepared submissions must not
be duplicated. Source identities and all eight cases are pinned in INPUTS.json;
the three branch plans declare exact geometry rules and denominators.

The shared adapter accepts these keys:

```text
branch: structure_informed_starts | solvent_guided | local_basin_breadth
inputs: {path, sha256} of common INPUTS.json
agreement: {path, sha256} of the branch plan
coordinate_limits: {maximum_angle_radian: 0.8, maximum_heavy_displacement_A: 0.8}
maximum_candidates_per_case: finite integer
cases: [{case_id, status, reason, candidates: [{id, full_q,
         optional coordinate: {path, sha256},
         optional native_reuse: {Ca: receipt_pin, La: receipt_pin}}]}]
```

The adapter checks the current common four-angle physical mapping, preserves
actual candidate XYZ coordinates after a 1e-12 Å map check, uses the existing
1e-7 Å final displacement tolerance, and deduplicates coordinates at 1e-12 Å.
It retains every old pool cell and cross-scores both metals for every admitted
new geometry. Branches retain unsupported proposal slots in their metadata;
failed required energy cells make the new pool unavailable, not baseline success.
The baseline always remains separately named.

Inspect the actual solvent-guided preparation without computing energies:

```bash
"$NIKASHA_PY" scripts/nikasha_finite_candidates.py validate \
  --manifest workspaces/solvent_guided_20260922/pool_v1/manifest.json
```

The existing `diagnostics/nikasha_shared_pool_20260922/run_mace.sbatch` and
`run_solvent.sbatch` execute finite manifests. Agent submission receipts record
the exact jobs, allocations and command arguments; consult them before any
execution. The source checkouts and historical protocols remain unchanged.

All three collections are finished. Replay their actual comparison without any
molecular evaluations. The temporary output directory is new on each invocation;
the archived result is never overwritten:

```bash
NIKASHA_REPLAY_DIR=$(mktemp -d /tmp/nikasha-parallel-replay.XXXXXX)
"$NIKASHA_PY" scripts/nikasha_parallel_compare.py \
  --inputs diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json \
  --collections \
    workspaces/structure_informed_starts_20260922/pool_v1/final_collection.json \
    workspaces/solvent_guided_20260922/pool_v1/after_solvent_0_1210101.json \
    workspaces/local_basin_breadth_20260922/scoring_v1/after_solvent_3_1210125.json \
  --output "$NIKASHA_REPLAY_DIR/COMPARISON.json"
```

The output preserves static and adaptive contrasts, individual Ca/La work,
selected candidates, and separate transfer of released/adaptive bands. It does
not calibrate a new threshold. The basin branch's conditional integral is
reported by its own analyzer; the shared comparison reports minima from the
expanded common pool (old candidates plus grid), not the standalone grid.
Its breadth term must not be added to an unrelated adaptive winner.

The primary archived result is
`workspaces/nikasha_parallel_pilots_20260922/COMPARISON_v1.json` and `.md`.
The separate integral is `workspaces/local_basin_breadth_20260922/RESULT_v2.json`.
Its earlier partialv1 is a preserved JSON serialization failure, not a scientific
output. Branch COMMANDS files give geometry, collection, integration and plot
replays. No scientific executable is needed for the comparison above.

All submitted jobs are terminal. Do not resubmit them from the generic runner
examples. Proposed native restarts, backend fallback, collective scaffold and
redox work have not run.

Validation logs: FINITE_TESTS.txt (5 passed), POOL_REGRESSION.txt (30 passed),
COMPARE_TESTS.txt (3 passed). These use actual pinned molecular fixtures and
explicit corrupted copies, and do not replace scientific integration runs.
