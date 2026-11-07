#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JAX_PLATFORMS=cpu
timeout --signal=TERM --kill-after=5s 180s srun -n1 -c3 --cpu-bind=cores bash -c '
workspaces/so3lr_tauri_screen_20261106/software/env/bin/python -u diagnostics/so3lr_neutral_pair_20261107/prepare_hydration.py &&
workspaces/so3lr_tauri_screen_20261106/software/env/bin/python -u diagnostics/so3lr_neutral_pair_20261107/trajectory.py --preparation workspaces/so3lr_neutral_pair_20261107/hydrated_preparation_v1/RESULT.json --output-dir "workspaces/so3lr_neutral_pair_20261107/hydrated_observations_${SLURM_JOB_ID}"'
