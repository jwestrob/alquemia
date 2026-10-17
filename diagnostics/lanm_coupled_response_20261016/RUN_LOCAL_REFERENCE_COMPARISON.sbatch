#!/bin/bash
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=lbt48_compare
set -euo pipefail
if (( $# != 3 )); then
    echo 'Usage: sbatch RUN_LOCAL_REFERENCE_COMPARISON.sbatch ABS_ADAPTER CONFIG FRESH_OUTPUT' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11 "$1" --config "$2" --output "$3"
