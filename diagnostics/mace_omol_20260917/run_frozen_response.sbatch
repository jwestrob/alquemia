#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-frozen-response
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=32G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_frozen_response.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=64
FROZEN_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FROZEN_ROOT=$(dirname "$1")
FROZEN_DRIVER=$FROZEN_ROOT/implementation/mace_frozen_response.py
"$FROZEN_PY" "$FROZEN_DRIVER" dry-run --manifest "$1"
/usr/bin/time -v -o "$FROZEN_ROOT/job_${SLURM_JOB_ID}.resources" \
    "$FROZEN_PY" "$FROZEN_DRIVER" execute --manifest "$1" --output "$FROZEN_ROOT/collection_job_$SLURM_JOB_ID.json"
