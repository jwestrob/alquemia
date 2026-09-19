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

## MACE proposal route (current)

Completed: proposal check1201847, all eight exact-rotation searches1201849.
Running: four DFT adjudication endpoints1201853. Preserve these tasks.
The MACE checks use actual first/third DFT gradient checkpoints and native OMOL;
they do not use masked features or assign free energies to water occupancies.

Replay preparation in fresh directories, without changing any source paths:

```bash
"$HYDRATION_PY" scripts/hydration_mace.py prepare \
  --manifests "$HYDRATION_ROOT/workspaces/hydration_network_20260918/orientation_1f6s_v1/manifest.json" \
              "$HYDRATION_ROOT/workspaces/hydration_network_20260918/orientation_6ip9_v1/manifest.json" \
  --software "$HYDRATION_ROOT/workspaces/mace_omol_20260917/software_v1/software_manifest.json" \
  --inventory "$HYDRATION_ROOT/workspaces/mace_canonical_20260916/audit_v2/inventory.json" \
  --agreement "$HYDRATION_ROOT/diagnostics/hydration_network_20260918/MACE_PROPOSAL_PLAN.md" \
  --output "$HYDRATION_ROOT/workspaces/hydration_network_20260918/mace_proposal_replay_v1"
"$HYDRATION_PY" scripts/hydration_proposal_opt.py prepare \
  --proposal-collection "$HYDRATION_ROOT/workspaces/hydration_network_20260918/mace_proposal_v2/collection_job_1201847.json" \
  --manifests "$HYDRATION_ROOT/workspaces/hydration_network_20260918/orientation_1f6s_v1/manifest.json" \
              "$HYDRATION_ROOT/workspaces/hydration_network_20260918/orientation_6ip9_v1/manifest.json" \
  --inventory "$HYDRATION_ROOT/workspaces/mace_canonical_20260916/audit_v2/inventory.json" \
  --agreement "$HYDRATION_ROOT/diagnostics/hydration_network_20260918/MACE_OPTIMIZATION_PLAN.md" \
  --output "$HYDRATION_ROOT/workspaces/hydration_network_20260918/mace_optimization_replay_v1"
```

For the optimization replay above, the existing finite executor command is:

```bash
sbatch --parsable --partition=gpu --nodelist=node-128-512g-8gpu-1 \
  --job-name=water-mace-opt --cpus-per-task=16 --gres=gpu:1 --mem=64474M \
  --export=ALL,MACE_MIN_MEMORY_MIB=64474 \
  --output="$HYDRATION_ROOT/workspaces/hydration_network_20260918/mace_optimization_replay_v1/slurm_%j.out" \
  --error="$HYDRATION_ROOT/workspaces/hydration_network_20260918/mace_optimization_replay_v1/slurm_%j.err" \
  diagnostics/mace_hybrid_20260916/run_pilot.sbatch \
  "$HYDRATION_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" \
  "$HYDRATION_ROOT/workspaces/hydration_network_20260918/mace_optimization_replay_v1/manifest.json" native
```

Do not rerun the completed searches without a scientific reason. Recollect
current results or prepare DFT replay from their validated output as follows:

```bash
"$HYDRATION_PY" scripts/hydration_proposal_opt.py collect \
  --manifest "$HYDRATION_ROOT/workspaces/hydration_network_20260918/mace_optimization_v1/manifest.json" \
  --output "$HYDRATION_ROOT/workspaces/hydration_network_20260918/mace_optimization_v1/recollection_v1.json"
"$HYDRATION_PY" scripts/hydration_adjudicate.py prepare \
  --collection "$HYDRATION_ROOT/workspaces/hydration_network_20260918/mace_optimization_v1/collection_job_1201849.json" \
  --agreement "$HYDRATION_ROOT/diagnostics/hydration_network_20260918/MACE_OPTIMIZATION_PLAN.md" \
  --output "$HYDRATION_ROOT/workspaces/hydration_network_20260918/proposal_dft_replay_v1"
```

DFT uses `adjudicate.sbatch` with the same absolute-manifest argument and log
options as `run.sbatch`. It runs four16-rank native energy/analytic-gradient
endpoints and writes its collection automatically. Actual proposal_dft_v1 had
one missing collector import copied from the pinned existing MACE implementation;
`collection_dependency_fix.json` records it. No quantum input or existing file
was changed. Future preparations include that dependency from the outset.

## Confirmed original-core result and archived GGR robustness

Completed job1201867, all four endpoints. The original recipe/SCF/core is retained;
only water H coordinates change. Read REPORT.md before applying the development
result elsewhere. Recollect and reproduce the expanded GGR comparison without
running any inference or DFT:

```bash
"$HYDRATION_PY" scripts/hydration_core_transfer.py collect \
  --manifest "$HYDRATION_ROOT/workspaces/hydration_network_20260918/core_transfer_v1/manifest.json" \
  --output "$HYDRATION_ROOT/workspaces/hydration_network_20260918/core_transfer_v1/recollection_v1.json"
"$HYDRATION_PY" scripts/hydration_core_transfer.py compare-ggr-replicates \
  --collection "$HYDRATION_ROOT/workspaces/hydration_network_20260918/core_transfer_v1/collection_1201867.json" \
  --ggr-study "$HYDRATION_ROOT/diagnostics/ggr_mechanism_plan_20260915/RESULT.json" \
  --output "$HYDRATION_ROOT/workspaces/hydration_network_20260918/core_transfer_v1/ggr_recomparison_v1.json"
```

For independently reproducing the preparation in a fresh directory:

```bash
"$HYDRATION_PY" scripts/hydration_core_transfer.py prepare \
  --config "$HYDRATION_ROOT/diagnostics/hydration_network_20260918/CONFIG.json" \
  --adjudication "$HYDRATION_ROOT/workspaces/hydration_network_20260918/proposal_dft_v1/collection_numeric_parser_v2.json" \
  --ggr-release "$HYDRATION_ROOT/diagnostics/baseline_benchmark_20260915/RESULT.json" \
  --agreement "$HYDRATION_ROOT/diagnostics/hydration_network_20260918/CORE_TRANSFER_PLAN.md" \
  --output "$HYDRATION_ROOT/workspaces/hydration_network_20260918/core_transfer_replay_v1"
```

Execution uses `core_transfer.sbatch` with the same explicit absolute-manifest
argument/log options as the other DFT wrappers. Do not submit a duplicate merely
to replay collection. The successful expanded DFT collection is
`proposal_dft_v1/collection_numeric_parser_v2.json`: its original collection
failed on the progress text `gCP correction ... done`, fixed by requiring actual
numeric tokens. All four quantum jobs were already successful; no recomputation.

The23 archived baseline/development tests can be rerun with:

```bash
"$HYDRATION_PY" -m unittest discover -s tests -p test_affordable_development.py -q
```

Current next action is to inspect1201824/1201825 and collect only when finished.
Do not infer final native convergence or occupancy from the confirmed preparation
result. The old DFT-only occupancy execution example above is historical prepared
capability, not the recommended continuation now that cheap MACE proposals work.
