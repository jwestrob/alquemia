#!/bin/bash
#SBATCH --job-name=Nikasha-Glu91-summary
#SBATCH --partition=standard-shared
#SBATCH --ntasks=1
#SBATCH --mem=0
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/glu91_response_20260928/COLLECT.py --config diagnostics/metal_environment_response_20260926/glu91_response_20260928/CONFIG.json --output workspaces/metal_environment_response_20260926/glu91_response_final_v1
