#!/bin/bash
#SBATCH --job-name=frozen-distribution
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=12
set -euo pipefail
if [ "$#" -lt 4 ]; then
    echo 'usage: RUN.sbatch canonical|transfer MANIFEST SHA OUTPUT [--calibration FILE --calibration-sha SHA]' >&2
    exit 2
fi
stage="$1"
manifest="$2"
sha="$3"
output="$4"
shift 4
runner="$(dirname "$manifest")/implementation/NATIVE.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
exec /usr/bin/python3 "$runner" "$stage" "$manifest" --sha "$sha" --output "$output" "$@"
