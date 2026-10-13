#!/bin/bash
#SBATCH --partition=gpu,gpu_h200
#SBATCH --job-name=Nikasha-Ca2-AF3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Slurm required}"
if [[ $# -ne 2 ]]; then echo 'Usage: FOLD.sbatch ABSOLUTE_MANIFEST_JSON ACCESSION' >&2; exit 2; fi
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK" MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export NUMEXPR_NUM_THREADS="$SLURM_CPUS_PER_TASK" VECLIB_MAXIMUM_THREADS="$SLURM_CPUS_PER_TASK" BLIS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export XLA_PYTHON_CLIENT_PREALLOCATE=false PYTHONUNBUFFERED=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/ca_coverage_recovery_20261012/refold/v3/RUN_AF3.py --manifest "$1" --target "$2"
