#!/bin/bash
#SBATCH -p gpu
#SBATCH -J triple100-search
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
TRIPLE_MANIFEST=$1
TRIPLE_SCRIPT=$2
TRIPLE_AGREEMENT=$3
for TRIPLE_PATH in "$TRIPLE_MANIFEST" "$TRIPLE_SCRIPT" "$TRIPLE_AGREEMENT"; do
 case "$TRIPLE_PATH" in /*) ;; *) exit 2;; esac
done
TRIPLE_RUN=$(dirname "$TRIPLE_MANIFEST")
TRIPLE_ROOT=$(dirname "$TRIPLE_RUN")
TRIPLE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$TRIPLE_RUN/implementation"
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" execute --manifest "$TRIPLE_MANIFEST"
TRIPLE_EXEC=$?
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" collect --manifest "$TRIPLE_MANIFEST" --output "$TRIPLE_RUN/collection.json" || exit $?
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" prepare_pool --proposals "$TRIPLE_RUN/collection.json" --agreement "$TRIPLE_AGREEMENT" --output "$TRIPLE_ROOT/pool" || exit $?
export PYTHONPATH="$TRIPLE_ROOT/pool/implementation"
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" execute_mace --manifest "$TRIPLE_ROOT/pool/manifest.json"
TRIPLE_CROSS=$?
if [ "$TRIPLE_EXEC" -ne 0 ]; then exit "$TRIPLE_EXEC"; fi
exit "$TRIPLE_CROSS"
