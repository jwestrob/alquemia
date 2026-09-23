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

After the three branch agents have finished their collections, generate a joint
comparison using their actual paths, for example the completed collection paths
recorded in CURRENT.md/REPORT.md. The command takes explicit inputs and refuses
to overwrite its output:

```bash
"$NIKASHA_PY" scripts/nikasha_parallel_compare.py --help
```

The output preserves static and adaptive contrasts, individual Ca/La work,
selected candidates, and separate transfer of released/adaptive bands. It does
not calibrate a new threshold. The basin branch's conditional integral is
reported by its own analyzer; the shared comparison only reports grid minima.

Validation logs: FINITE_TESTS.txt (5 passed), POOL_REGRESSION.txt (30 passed),
COMPARE_TESTS.txt (3 passed). These use actual pinned molecular fixtures and
explicit corrupted copies, and do not replace scientific integration runs.
