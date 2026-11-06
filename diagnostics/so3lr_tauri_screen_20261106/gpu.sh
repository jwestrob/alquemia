#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false JAX_ENABLE_X64=true
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
workspaces/so3lr_tauri_screen_20261106/software/env/bin/python -u diagnostics/so3lr_tauri_screen_20261106/run_tauri.py "$@"
