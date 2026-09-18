#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-GK-source-check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64000M
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_gk_source_native.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=64
SOURCE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
SOURCE_ROOT=$(dirname "$1")
SOURCE_DRIVER=$SOURCE_ROOT/implementation/mace_gk_source_native.py
"$SOURCE_PY" "$SOURCE_DRIVER" dry-run --manifest "$1"
"$SOURCE_PY" "$SOURCE_DRIVER" execute --manifest "$1"
"$SOURCE_PY" "$SOURCE_DRIVER" collect --manifest "$1" --output "$SOURCE_ROOT/collection_job_$SLURM_JOB_ID.json"
