#!/bin/bash
#SBATCH -p gpu
#SBATCH -J envelope100-search
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -euo pipefail
ENVELOPE_SHARD=$1
case "$ENVELOPE_SHARD" in /*) ;; *) exit 2;; esac
ENVELOPE_CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
ENVELOPE_IMPL="$ENVELOPE_SHARD/origins/implementation"
ENVELOPE_AGREEMENT=$2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$ENVELOPE_CPU" "$ENVELOPE_IMPL/motion_envelope_transfer.py" collect_origins --stage "$ENVELOPE_SHARD/origins" --output "$ENVELOPE_SHARD/origins/COLLECTION.json"
"$ENVELOPE_CPU" "$ENVELOPE_IMPL/motion_envelope_transfer.py" prepare_searches --origins "$ENVELOPE_SHARD/origins/COLLECTION.json" --output "$ENVELOPE_SHARD/searches"
set +e
"$ENVELOPE_CPU" "$ENVELOPE_SHARD/searches/implementation/motion_envelope_transfer.py" execute_searches --manifest "$ENVELOPE_SHARD/searches/manifest.json"
ENVELOPE_RC=$?
set -e
"$ENVELOPE_CPU" "$ENVELOPE_SHARD/searches/implementation/motion_envelope_transfer.py" collect_searches --manifest "$ENVELOPE_SHARD/searches/manifest.json" --output "$ENVELOPE_SHARD/searches/COLLECTION.json"
"$ENVELOPE_CPU" "$ENVELOPE_SHARD/searches/implementation/motion_envelope_transfer.py" prepare_pool --proposals "$ENVELOPE_SHARD/searches/COLLECTION.json" --agreement "$ENVELOPE_AGREEMENT" --output "$ENVELOPE_SHARD/pool"
exit "$ENVELOPE_RC"
