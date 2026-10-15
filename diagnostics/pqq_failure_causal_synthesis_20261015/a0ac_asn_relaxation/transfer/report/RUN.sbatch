#!/bin/bash
#SBATCH --job-name=Asn-transfer-report
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /usr/bin/python3 "$1" --collection "$2" --execution "$3" --output "$4"
