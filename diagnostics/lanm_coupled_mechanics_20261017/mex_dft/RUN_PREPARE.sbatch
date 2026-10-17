#!/bin/bash
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --nodelist=node-64-768g-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=mex_dft_prepare
set -euo pipefail
if (( $# != 4 )); then
    echo 'Usage: sbatch RUN_PREPARE.sbatch ABS_ADAPTER CONFIG FRESH_OUTPUT prepare|collect' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export BLIS_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE PYTHONUNBUFFERED=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11 "$1" --config "$2" --output "$3" --stage "$4"
