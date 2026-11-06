#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout --signal=TERM --kill-after=2s 30s \
  workspaces/so3lr_tauri_screen_20261106/software/env/bin/python \
  diagnostics/so3lr_swap_spectator_20261106/analyze.py "$@"
