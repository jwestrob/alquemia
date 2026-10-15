#!/bin/bash
#SBATCH --job-name=MMOL-common3-comparison
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue
set -euo pipefail
REFERENCE_MANIFEST="$1"
REFERENCE_SHA="$2"
REFERENCE_OUTPUT="$3"
shift 3
REFERENCE_WORK=$(dirname "$REFERENCE_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE LC_ALL=C LANG=C
exec /usr/bin/python3 "$REFERENCE_WORK/implementation/COMPARE.py" "$REFERENCE_MANIFEST" --manifest-sha "$REFERENCE_SHA" --output "$REFERENCE_OUTPUT" --collections "$@"
