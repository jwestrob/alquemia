#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-distributed-source
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_distributed_source.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
SOURCE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
SOURCE_ROOT=$(dirname "$1")
SOURCE_DRIVER=$SOURCE_ROOT/implementation/mace_distributed_source.py
"$SOURCE_PY" "$SOURCE_DRIVER" execute --manifest "$1" \
    --output "$SOURCE_ROOT/execution_job_${SLURM_JOB_ID}.json"
"$SOURCE_PY" "$SOURCE_DRIVER" collect --manifest "$1" \
    --output "$SOURCE_ROOT/report_job_${SLURM_JOB_ID}"
