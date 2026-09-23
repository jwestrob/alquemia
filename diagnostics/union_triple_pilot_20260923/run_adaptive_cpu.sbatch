#!/bin/bash
#SBATCH -p gpu
#SBATCH -J triple-pool-gfn
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=65536M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
TRIPLE_ROOT=$1
case "$TRIPLE_ROOT" in /*) ;; *) exit 2;; esac
TRIPLE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
TRIPLE_POOL="$TRIPLE_ROOT/pool_v1"
TRIPLE_SCRIPT="$TRIPLE_POOL/implementation/union_triple_adaptive.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" validate_pool --manifest "$TRIPLE_POOL/manifest.json" || exit $?
"$TRIPLE_CPU" "$TRIPLE_POOL/implementation/affordable_workflow.py" execute --manifest "$TRIPLE_POOL/solvent/shard_0/manifest.json"
TRIPLE_EXEC=$?
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" collect_pool --manifest "$TRIPLE_POOL/manifest.json" --output "$TRIPLE_POOL/collection_final.json" || exit $?
"$TRIPLE_CPU" "$TRIPLE_SCRIPT" reference --collection "$TRIPLE_POOL/collection_final.json" --output "$TRIPLE_ROOT/REFERENCE_v1.json"
TRIPLE_REPORT=$?
if [ "$TRIPLE_EXEC" -ne 0 ]; then exit "$TRIPLE_EXEC"; fi
exit "$TRIPLE_REPORT"
