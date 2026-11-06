#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JAX_PLATFORMS=cpu JAX_ENABLE_X64=true
srun -n 1 -c 12 --cpu-bind=cores workspaces/so3lr_tauri_screen_20261106/software/env/bin/python -u diagnostics/so3lr_tauri_screen_20261106/run_tauri.py --cpu "$@"
