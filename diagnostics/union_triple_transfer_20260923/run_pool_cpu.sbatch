#!/bin/bash
#SBATCH -p gpu,standard-shared
#SBATCH -J triple100-pool
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=65536M
set -uo pipefail
TRIPLE_ROOT=$1
TRIPLE_SCRIPT=$2
case "$TRIPLE_ROOT" in /*) ;; *) exit 2;; esac
case "$TRIPLE_SCRIPT" in /*) ;; *) exit 2;; esac
TRIPLE_POOL="$TRIPLE_ROOT/pool"
TRIPLE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONPATH="$TRIPLE_POOL/implementation"
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" validate_pool --manifest "$TRIPLE_POOL/manifest.json" || exit $?
"$TRIPLE_CPU" "$TRIPLE_POOL/implementation/affordable_workflow.py" execute --manifest "$TRIPLE_POOL/solvent/shard_0/manifest.json"
TRIPLE_EXEC=$?
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" collect_pool --manifest "$TRIPLE_POOL/manifest.json" --output "$TRIPLE_POOL/collection_final.json" || exit $?
exit "$TRIPLE_EXEC"
