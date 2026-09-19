# Joint occupancy operations

From the repository root:

```bash
HYDRATION_ROOT=$PWD
HYDRATION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

The source inventory is prepared at
`workspaces/hydration_occupancy_20260918/states_v1/manifest.json`.
Its36 native-optimization input templates are **preparation inputs only**;
never submit that source manifest. The MACE stage consumes the same geometry
and prepares32 short rigid-water searches in16 finite executor tasks.
Actual submitted MACE job:1201908. Do not duplicate completed/running tasks.

## Preparation replay, in fresh directories

```bash
"$HYDRATION_PY" scripts/hydration_network.py prepare \
  --config diagnostics/hydration_network_20260918/CONFIG.json \
  --output workspaces/hydration_occupancy_20260918/states_replay_v1 \
  --occupancy --exclude-full
"$HYDRATION_PY" scripts/hydration_proposal_opt.py prepare --occupancy \
  --proposal-collection workspaces/hydration_network_20260918/mace_proposal_v2/collection_job_1201847.json \
  --manifests workspaces/hydration_occupancy_20260918/states_replay_v1/manifest.json \
  --inventory workspaces/mace_canonical_20260916/audit_v2/inventory.json \
  --agreement diagnostics/hydration_occupancy_20260918/AGREEMENT.md \
  --output workspaces/hydration_occupancy_20260918/mace_replay_v1
```

The existing GPU executor is used with one16CPU/oneGPU/64474MiB allocation:

```bash
sbatch --parsable --partition=gpu --nodelist=node-128-512g-8gpu-1 \
  --job-name=water-occupancy-mace --cpus-per-task=16 --gres=gpu:1 --mem=64474M \
  --export=ALL,MACE_MIN_MEMORY_MIB=64474 \
  --output="$HYDRATION_ROOT/workspaces/hydration_occupancy_20260918/mace_replay_v1/slurm_%j.out" \
  --error="$HYDRATION_ROOT/workspaces/hydration_occupancy_20260918/mace_replay_v1/slurm_%j.err" \
  diagnostics/mace_hybrid_20260916/run_pilot.sbatch \
  "$HYDRATION_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python" \
  "$HYDRATION_ROOT/workspaces/hydration_occupancy_20260918/mace_replay_v1/manifest.json" native
```

## Actual DFT continuation

After job1201908 is complete and all16 proposals qualify:

```bash
"$HYDRATION_PY" scripts/hydration_occupancy.py prepare \
  --states workspaces/hydration_occupancy_20260918/states_v1/manifest.json \
  --proposals workspaces/hydration_occupancy_20260918/mace_v1/collection_job_1201908.json \
  --full-collection workspaces/hydration_network_20260918/proposal_dft_v1/collection_numeric_parser_v2.json \
  --agreement diagnostics/hydration_occupancy_20260918/AGREEMENT.md \
  --output workspaces/hydration_occupancy_20260918/dft_v1
"$HYDRATION_PY" scripts/affordable_workflow.py dry-run \
  --manifest workspaces/hydration_occupancy_20260918/dft_v1/manifest.json
sbatch --parsable \
  --output="$HYDRATION_ROOT/workspaces/hydration_occupancy_20260918/dft_v1/slurm_%j.out" \
  --error="$HYDRATION_ROOT/workspaces/hydration_occupancy_20260918/dft_v1/slurm_%j.err" \
  diagnostics/hydration_occupancy_20260918/run_dft.sbatch \
  "$HYDRATION_ROOT/workspaces/hydration_occupancy_20260918/dft_v1/manifest.json"
```

This runs20 new native endpoints (four16-rank workers); four full states reuse
actual matching results. Source graph/coordinates, charges, method, state coverage
and receipts are checked explicitly. No MACE energies compare different water
counts. Completed output and failed attempts are retained by the existing runner.

## Collection and comparison

After execution, the wrapper writes its collection automatically. A fresh manual
recollection and reference comparison can be made with:

```bash
"$HYDRATION_PY" scripts/hydration_occupancy.py collect \
  --manifest workspaces/hydration_occupancy_20260918/dft_v1/manifest.json \
  --output workspaces/hydration_occupancy_20260918/dft_v1/recollection_v1.json
"$HYDRATION_PY" scripts/hydration_occupancy.py analyze \
  --collection workspaces/hydration_occupancy_20260918/dft_v1/recollection_v1.json \
  --reference workspaces/hydration_network_20260918/water_reference_v2/reference_1201831.json \
  --output diagnostics/hydration_occupancy_20260918/RESULT_recollection_v1.json
```

Paths are write-once. Use a fresh explicit output path when repeating collection.
Incomplete states remain unavailable; no bound free-energy correction, occupancy
probability or calibrated decision is generated from bare electronic energies.

```bash
"$HYDRATION_PY" -m unittest discover -s tests -p test_hydration_occupancy.py -v
"$HYDRATION_PY" -m unittest discover -s tests -p test_hydration_network.py -v
```

## Physical-gradient diagnostic and report exports

These operations use saved scientific outputs and launch no new calculations:

```bash
"$HYDRATION_PY" scripts/hydration_occupancy.py motion \
  --collection workspaces/hydration_occupancy_20260918/dft_v1/recollection_v1.json \
  --output diagnostics/hydration_occupancy_20260918/MOTION_recollection_v1.json
"$HYDRATION_PY" scripts/hydration_occupancy.py report \
  --result diagnostics/hydration_occupancy_20260918/RESULT_recollection_v1.json \
  --output diagnostics/hydration_occupancy_20260918/export_recollection_v1
"$HYDRATION_PY" scripts/hydration_site_proposals.py \
  --states workspaces/hydration_occupancy_20260918/states_v1/manifest.json \
  --agreement diagnostics/hydration_occupancy_20260918/SITE_PROPOSAL_PLAN.md \
  --output diagnostics/hydration_occupancy_20260918/SITE_PROPOSALS_recollection_v1.json
```

`report` exports CSV, Markdown, SVG and PNG. The projected gradients are physical
water translations/rotations, not thermal fluctuations or a validated curvature.
Experimental-site proposals must also be checked for water-water incompatibility:
the unmatched transferred site in this pair is1.912A from an existing water and
cannot simply be co-occupied at its proposed position. No site was inserted.
