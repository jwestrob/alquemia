#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=nikasha-mxaf-recovery
#SBATCH --nodes=1
#SBATCH --ntasks=7
#SBATCH --no-requeue
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/ca_coverage_recovery_20261012/mxaf_v2/slurm_%j.out
set -euo pipefail
test "$#" -eq 1
test -n "${SLURM_JOB_ID:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -u diagnostics/ca_coverage_recovery_20261012/mxaf/RECOVER.py execute --manifest "$1"
