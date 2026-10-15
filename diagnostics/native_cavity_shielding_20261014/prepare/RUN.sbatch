#!/bin/bash
#SBATCH --job-name=native-cavity-prepare
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 "$(dirname "$1")/implementation/PREPARE.py" execute --manifest "$1" --manifest-sha "$2"
