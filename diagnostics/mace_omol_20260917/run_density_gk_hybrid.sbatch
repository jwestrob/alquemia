#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-density-gk-hybrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64000M
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_density_gk_hybrid.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=64
DENSITY_GK_ROOT=$(dirname "$1")
DENSITY_GK_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
DENSITY_GK_DRIVER=$DENSITY_GK_ROOT/implementation/mace_density_gk_hybrid.py
"$DENSITY_GK_PY" "$DENSITY_GK_DRIVER" dry-run --manifest "$1"
"$DENSITY_GK_PY" "$DENSITY_GK_DRIVER" execute --manifest "$1"
"$DENSITY_GK_PY" "$DENSITY_GK_DRIVER" collect --manifest "$1" \
    --output "$DENSITY_GK_ROOT/collection_job_${SLURM_JOB_ID}.json"
