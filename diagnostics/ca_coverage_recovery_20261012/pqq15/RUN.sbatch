#!/bin/bash
#SBATCH --partition=standard
#SBATCH --job-name=Nikasha-PQQ15-recovery
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Preparation requires Slurm}"
if [[ $# -ne 1 ]]; then echo 'Usage: RUN.sbatch ABSOLUTE_MANIFEST_JSON' >&2; exit 2; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/ca_coverage_recovery_20261012/pqq15/PREPARE.py execute --manifest "$1"
