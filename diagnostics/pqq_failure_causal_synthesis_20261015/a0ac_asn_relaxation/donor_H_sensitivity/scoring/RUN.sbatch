#!/bin/bash
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=7
#SBATCH --job-name=Asn-twoH-MACE
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual Slurm required}"
test "$#" -eq 3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python "$(dirname "$1")/SCORE.py" "$3" --config "$1" --config-sha "$2"
