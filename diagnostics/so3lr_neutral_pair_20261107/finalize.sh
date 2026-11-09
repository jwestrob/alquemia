#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 JAX_PLATFORMS=cpu
/groups/banfield/users/jwestrob/envs/lnsel/bin/python diagnostics/so3lr_neutral_pair_20261107/finalize.py \
 --output-dir "workspaces/so3lr_neutral_pair_20261107/final_review_${SLURM_JOB_ID}"
