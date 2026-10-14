#!/bin/bash
#SBATCH --job-name=frozen-distribution-score
#SBATCH --partition=standard
#SBATCH --ntasks=1
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ "$#" -lt 1 ]; then
    echo 'usage: RUN.sbatch FROZEN_SCORE.py --config CONFIG --sha SHA --stage canonical|transfer --native-collection FILE [--canonical-calibration FILE] --output FILE' >&2
    exit 2
fi
exec /usr/bin/python3 "$@"
