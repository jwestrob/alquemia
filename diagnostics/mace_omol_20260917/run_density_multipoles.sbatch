#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-density-multipoles
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_density_multipoles.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
MULTIPOLE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
MULTIPOLE_ROOT=$(dirname "$1")
MULTIPOLE_DRIVER=$MULTIPOLE_ROOT/implementation/mace_density_multipoles.py
"$MULTIPOLE_PY" "$MULTIPOLE_DRIVER" execute --manifest "$1" \
    --output "$MULTIPOLE_ROOT/execution_job_${SLURM_JOB_ID}.json"
"$MULTIPOLE_PY" "$MULTIPOLE_DRIVER" collect --manifest "$1" \
    --output "$MULTIPOLE_ROOT/report_job_${SLURM_JOB_ID}"
