#!/bin/bash
# Saved native trajectory only; no potential or simulation calls.
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JAX_PLATFORMS=cpu
workspaces/so3lr_tauri_screen_20261106/software/env/bin/python diagnostics/so3lr_neutral_pair_20261107/trajectory.py \
 --md-dir workspaces/so3lr_neutral_pair_20261107/md_scout_v1 \
 --preparation workspaces/so3lr_neutral_pair_20261107/hydrated_preparation_v1/RESULT.json \
 --output-dir "workspaces/so3lr_neutral_pair_20261107/md_scout_analysis_${SLURM_JOB_ID}"
