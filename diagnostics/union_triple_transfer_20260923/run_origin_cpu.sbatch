#!/bin/bash
#SBATCH -p gpu,standard-shared
#SBATCH -J triple100-q0-gfn
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=65536M
set -uo pipefail
TRIPLE_STAGE=$1
TRIPLE_SCRIPT=$2
case "$TRIPLE_STAGE" in /*) ;; *) exit 2;; esac
case "$TRIPLE_SCRIPT" in /*) ;; *) exit 2;; esac
TRIPLE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$TRIPLE_STAGE/implementation"
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" execute_origins --manifest "$TRIPLE_STAGE/solvent/manifest.json"
TRIPLE_EXEC=$?
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" collect_origins --stage "$TRIPLE_STAGE" --output "$TRIPLE_STAGE/COLLECTION.json" || exit $?
exit "$TRIPLE_EXEC"
