#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
task_work="$PWD/workspaces/so3lr_neutral_pair_20261107"
timeout --signal=TERM --kill-after=3s 180s make -C "$task_work/tools/packmol_src" -j1 > "$task_work/packmol_build.log" 2>&1
timeout --signal=TERM --kill-after=5s 600s /groups/banfield/users/jwestrob/envs/lnsel/bin/python -u diagnostics/so3lr_neutral_pair_20261107/prepare.py
