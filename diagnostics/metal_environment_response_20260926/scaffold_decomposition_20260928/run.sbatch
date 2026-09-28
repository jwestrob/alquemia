#!/bin/bash
#SBATCH --job-name=Nikasha-scaffold-account
#SBATCH --partition=standard-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --no-requeue
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ "$(date +%s)" -ge "$(TZ=America/Los_Angeles date -d '2026-09-28 10:00:00' +%s)" ]; then exit 75; fi
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/scaffold_decomposition_20260928/RUN.py workspaces/metal_environment_response_20260926/scaffold_decomposition_v1
