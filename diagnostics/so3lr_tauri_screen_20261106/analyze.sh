#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT=$(git rev-parse --show-toplevel)
exec "$ROOT/workspaces/so3lr_tauri_screen_20261106/software/env/bin/python" "$ROOT/diagnostics/so3lr_tauri_screen_20261106/analyze_tauri.py" "$@"
