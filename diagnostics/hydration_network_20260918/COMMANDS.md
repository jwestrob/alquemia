# Hydration-network operations

Run from this repository. These are development operations; production is unchanged.
Set explicit paths once:

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
HYDRATION_ROOT=$PWD
HYDRATION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Current work: do not resubmit

Jobs 1201824/1201825 run the eight full-water orientation searches. The shared
water reference is complete as1201831. Recover [collection policy](COLLECTION_POLICY.md):
the protein jobs' original runner expects an unnecessary `.engrad`, so use the
corrected collector after they terminate, even if Slurm labels the jobs failed
for that artifact alone. Nonconverged quantum optimizations remain unavailable.

```bash
squeue -j 1201824,1201825
"$HYDRATION_PY" scripts/hydration_network.py collect \
  --manifest "$HYDRATION_ROOT/workspaces/hydration_network_20260918/orientation_1f6s_v1/manifest.json" \
  --output "$HYDRATION_ROOT/workspaces/hydration_network_20260918/orientation_1f6s_v1/collected_opt_v2.json"
"$HYDRATION_PY" scripts/hydration_network.py collect \
  --manifest "$HYDRATION_ROOT/workspaces/hydration_network_20260918/orientation_6ip9_v1/manifest.json" \
  --output "$HYDRATION_ROOT/workspaces/hydration_network_20260918/orientation_6ip9_v1/collected_opt_v2.json"
"$HYDRATION_PY" scripts/hydration_state_analysis.py \
  --collections "$HYDRATION_ROOT/workspaces/hydration_network_20260918/orientation_1f6s_v1/collected_opt_v2.json" \
                "$HYDRATION_ROOT/workspaces/hydration_network_20260918/orientation_6ip9_v1/collected_opt_v2.json" \
  --reference "$HYDRATION_ROOT/workspaces/hydration_network_20260918/water_reference_v2/reference_1201831.json" \
  --previous "$HYDRATION_ROOT/workspaces/hydration_square_20260918/repaired_v1/collection_1201801.json" \
  --output "$HYDRATION_ROOT/diagnostics/hydration_network_20260918/ORIENTATION_RESULT.json"
```

Output paths are write-once. Do not collect prematurely into these final paths;
use a separately named checkpoint if inspecting an incomplete job.

## Reproducible preparation / future occupancy continuation

The preparer audits pinned parents, discovers source-defined neighbors, completes
chemical fragments and verifies paired coordinates/charges. `--case` optionally
selects 1F6S or6IP9 while retaining the common protein fragment selection.
The unsubmitted `occupancy_prepared_v1` has the obsolete artifact expectation.
After assessing Stage A, prepare its replacement without duplicating full states:

```bash
"$HYDRATION_PY" scripts/hydration_network.py prepare \
  --config "$HYDRATION_ROOT/diagnostics/hydration_network_20260918/CONFIG.json" \
  --output "$HYDRATION_ROOT/workspaces/hydration_network_20260918/occupancy_v2" \
  --occupancy --exclude-full
"$HYDRATION_PY" scripts/affordable_workflow.py dry-run \
  --manifest "$HYDRATION_ROOT/workspaces/hydration_network_20260918/occupancy_v2/manifest.json"
```

This prepares36 tasks (32 constrained orientation optimizations, four empty-state
single points). Preparation is not submission. The existing finite-manifest
executor is used; no allocated compute-time stopping rule is applied:

```bash
sbatch \
  --output="$HYDRATION_ROOT/workspaces/hydration_network_20260918/occupancy_v2/slurm_%j.out" \
  --error="$HYDRATION_ROOT/workspaces/hydration_network_20260918/occupancy_v2/slurm_%j.err" \
  diagnostics/hydration_network_20260918/run.sbatch \
  "$HYDRATION_ROOT/workspaces/hydration_network_20260918/occupancy_v2/manifest.json"
```

The runner uses four16-rank workers on64CPUs. Collection is written automatically
by this corrected preparation. Combine that collection with the two Stage A
collections using the same analysis command, with a fresh result path. Do not
use an incomplete state table to choose occupancies or a favorable deletion.

## Tests

```bash
"$HYDRATION_PY" -m unittest discover -s tests -p test_hydration_network.py -v
"$HYDRATION_PY" -m unittest discover -s tests -p test_hydration_square.py -v
```

Fixture absence causes explicit skips; the completed water reference is an actual
scientific integration result, not a mock. Missing bound-state entropy and
non-electrostatic terms remain unavailable. No new classification bands.
