#!/bin/bash
#SBATCH --job-name=canonical-saved-summary
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
REPORT_DIR=$(dirname "$1")
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11 "$REPORT_DIR/implementation/SUMMARY.py" execute "$1" "$2"
