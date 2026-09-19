# Local water-motion operations

Run from the repository root. All outputs are write-once. Use fresh explicit
paths for replay; do not resubmit completed jobs.

```bash
HYDRATION_ROOT=$PWD
HYDRATION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Prepare and validate

This prepares exactly eight DFT and sixteen MACE evaluations on four recorded
centers, with no optimization or change in water count. The implemented science
and acceptance limits are in [AGREEMENT.md](AGREEMENT.md).

```bash
"$HYDRATION_PY" scripts/hydration_water_motion.py prepare \
  --collection workspaces/hydration_occupancy_20260918/dft_v1/collection_1201910.json \
  --inventory workspaces/mace_canonical_20260916/audit_v2/inventory.json \
  --agreement diagnostics/hydration_motion_20260918/AGREEMENT.md \
  --output workspaces/hydration_motion_20260918/replay_v1
"$HYDRATION_PY" scripts/affordable_workflow.py dry-run \
  --manifest workspaces/hydration_motion_20260918/replay_v1/dft/manifest.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  workspaces/mace_hybrid_20260916/software_v1/venv/bin/python \
  workspaces/hydration_motion_20260918/replay_v1/mace/implementation/mace_hybrid.py dry-run \
  --manifest workspaces/hydration_motion_20260918/replay_v1/mace/manifest.json
```

## Execute through existing runners

These are replay commands, not an instruction to rerun completed calculations.
The actual DFT job1201953 uses `pilot_v2/dft`; actual MACE job1201958 uses
`mace_retry_v1`. Failed MACE startup1201954 ran zero inference calls; its1s
allocation remains included in cost. Both successful manifests share identical
physical preparation. Preparation-only `pilot_v1` was never submitted.

```bash
sbatch --parsable \
  --output="$HYDRATION_ROOT/workspaces/hydration_motion_20260918/replay_v1/dft/slurm_%j.out" \
  --error="$HYDRATION_ROOT/workspaces/hydration_motion_20260918/replay_v1/dft/slurm_%j.err" \
  diagnostics/hydration_motion_20260918/run_dft.sbatch \
  "$HYDRATION_ROOT/workspaces/hydration_motion_20260918/replay_v1/dft/manifest.json"
sbatch --parsable --partition=gpu --nodelist=node-128-512g-8gpu-1 \
  --job-name=water-motion-mace --cpus-per-task=16 --gres=gpu:1 --mem=64474M \
  --export=ALL,MACE_MIN_MEMORY_MIB=64474 \
  --output="$HYDRATION_ROOT/workspaces/hydration_motion_20260918/replay_v1/mace/slurm_%j.out" \
  --error="$HYDRATION_ROOT/workspaces/hydration_motion_20260918/replay_v1/mace/slurm_%j.err" \
  diagnostics/mace_hybrid_20260916/run_pilot.sbatch \
  "$HYDRATION_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" \
  "$HYDRATION_ROOT/workspaces/hydration_motion_20260918/replay_v1/mace/manifest.json" native
```

## Collect and compare actual results

Wrappers collect automatically. This read-only scientific replay re-parses actual
endpoint outputs and gradients, validates cached artifacts, and compares the
DFT-anchored cheap energy/gradient changes. It submits no calculation.

```bash
"$HYDRATION_PY" scripts/hydration_water_motion.py collect-dft \
  --manifest workspaces/hydration_motion_20260918/pilot_v2/dft/manifest.json \
  --output workspaces/hydration_motion_20260918/pilot_v2/dft/recollection_v1.json
"$HYDRATION_PY" scripts/mace_hybrid.py collect \
  --manifest workspaces/hydration_motion_20260918/mace_retry_v1/manifest.json \
  --output workspaces/hydration_motion_20260918/mace_retry_v1/recollection_v1.json
"$HYDRATION_PY" scripts/hydration_water_motion.py compare \
  --preparation workspaces/hydration_motion_20260918/pilot_v2/preparation.json \
  --mace-collection workspaces/hydration_motion_20260918/mace_retry_v1/recollection_v1.json \
  --dft-collection workspaces/hydration_motion_20260918/pilot_v2/dft/recollection_v1.json \
  --output diagnostics/hydration_motion_20260918/RESULT_recollection_v1.json
```

`refresh-mace --manifest PATH --output FRESH_DIRECTORY` creates a new executor
snapshot for a technical retry, preserving the same geometry/preparation/model
and the previous manifest. No scientific result is silently overwritten or reused
across changed implementations. Run its pinned dry-run with the GPU environment
before submission. Failed preparation/allocations stay visible.

```bash
"$HYDRATION_PY" -m unittest discover -s tests -p 'test_hydration*.py' -v
"$HYDRATION_PY" -m unittest discover -s tests -p test_affordable_development.py -v
```

No occupancy probabilities, entropy or relaxation score is released by this
operation: a single physical coordinate cannot validate a complete water basin.

## Export the measured result

```bash
"$HYDRATION_PY" scripts/hydration_water_motion.py report \
  --result diagnostics/hydration_motion_20260918/RESULT_v2.json \
  --output diagnostics/hydration_motion_20260918/export_replay_v1
```

The CSV/table retain all four endpoints. The energy plot uses actual sampled
energies and the declared anchored prediction; connecting lines guide the eye
and are not a fitted energy surface.

The final numerical record is `RESULT_v2.json`, reproduced with the frozen
`workspaces/hydration_motion_20260918/analysis_implementation_v1/` modules.
Its scientific contents equal the first `RESULT.json`; only the analyzer source
pin now points to an immutable snapshot. Final exports are `export_v3/`.
Earlier plots are retained; v1 omitted the reused central anchor from the
connecting line, corrected in v2/v3 without changing any energy.
