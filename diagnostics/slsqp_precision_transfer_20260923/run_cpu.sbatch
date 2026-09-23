#!/bin/bash
#SBATCH -p gpu
#SBATCH -J precision225-gfn
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
TRANSFER_MANIFEST=$1
case "$TRANSFER_MANIFEST" in /*) ;; *) exit 2;; esac
TRANSFER_POOL=$(dirname "$TRANSFER_MANIFEST")
TRANSFER_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$TRANSFER_PY" "$TRANSFER_POOL/implementation/affordable_workflow.py" execute --manifest "$TRANSFER_POOL/solvent/shard_0/manifest.json"
TRANSFER_EXEC=$?
"$TRANSFER_PY" "$TRANSFER_POOL/implementation/slsqp_precision_transfer.py" collect_pool --manifest "$TRANSFER_MANIFEST" --output "$TRANSFER_POOL/collection_final.json"
TRANSFER_COLLECT=$?
if [ "$TRANSFER_EXEC" -ne 0 ]; then exit "$TRANSFER_EXEC"; fi
exit "$TRANSFER_COLLECT"
