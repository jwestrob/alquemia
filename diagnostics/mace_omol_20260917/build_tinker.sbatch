#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-tinker-build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64000M
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: build_tinker.sbatch ABSOLUTE_BUILD_MANIFEST' >&2
    exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
TINKER_BUILD_ROOT=$(dirname "$1")
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
    "$TINKER_BUILD_ROOT/build.py" "$1"
