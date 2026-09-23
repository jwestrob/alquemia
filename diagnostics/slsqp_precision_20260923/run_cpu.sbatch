#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-ftol-gfn
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --nodelist=node-128-512g-8gpu-1
set -uo pipefail
PRECISION_MANIFEST=$1
case "$PRECISION_MANIFEST" in /*) ;; *) exit 2;; esac
PRECISION_POOL=$(dirname "$PRECISION_MANIFEST")
PRECISION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$PRECISION_PY" "$PRECISION_POOL/implementation/affordable_workflow.py" execute --manifest "$PRECISION_POOL/solvent/shard_0/manifest.json"
PRECISION_EXEC=$?
"$PRECISION_PY" "$PRECISION_POOL/implementation/slsqp_precision.py" collect_pool --manifest "$PRECISION_MANIFEST" --output "$PRECISION_POOL/collection_final.json"
PRECISION_COLLECT=$?
if [ "$PRECISION_EXEC" -ne 0 ]; then exit "$PRECISION_EXEC"; fi
exit "$PRECISION_COLLECT"
