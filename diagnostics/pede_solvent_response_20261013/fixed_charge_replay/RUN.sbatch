#!/bin/bash
#SBATCH --job-name=PedE-fixed-charge-ALPB
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --no-requeue
set -euo pipefail
if [ "$#" -lt 4 ]; then echo 'Usage: RUN.sbatch build-test-execute|build|test|execute ABS_MANIFEST SHA NEW_OUTPUT [reuse options]' >&2; exit 2; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
replay_action=$1
replay_manifest=$2
replay_sha=$3
replay_output=$4
shift 4
exec /usr/bin/python3 "$(dirname "$replay_manifest")/implementation/REPLAY.py" "$replay_action" "$replay_manifest" --sha "$replay_sha" --output "$replay_output" "$@"
