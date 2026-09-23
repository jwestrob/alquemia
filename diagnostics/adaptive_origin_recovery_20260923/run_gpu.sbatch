#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-origin-adaptive
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
RECOVERY_MANIFEST=$1
case "$RECOVERY_MANIFEST" in /*) ;; *) exit 2;; esac
RECOVERY_RUN=$(dirname "$RECOVERY_MANIFEST")
RECOVERY_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RECOVERY_SCRIPT="$RECOVERY_RUN/implementation/adaptive_origin_recovery.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$RECOVERY_PY" "$RECOVERY_SCRIPT" execute --manifest "$RECOVERY_MANIFEST"
RECOVERY_EXEC=$?
"$RECOVERY_PY" "$RECOVERY_SCRIPT" prepare-pool --proposals "$RECOVERY_MANIFEST" --output "$RECOVERY_RUN/pool"
RECOVERY_PREP=$?
if [ "$RECOVERY_PREP" -ne 0 ]; then exit "$RECOVERY_PREP"; fi
"$RECOVERY_PY" "$RECOVERY_RUN/pool/implementation/adaptive_origin_recovery.py" execute-cross --manifest "$RECOVERY_RUN/pool/manifest.json"
RECOVERY_CROSS=$?
if [ "$RECOVERY_EXEC" -ne 0 ]; then exit "$RECOVERY_EXEC"; fi
exit "$RECOVERY_CROSS"
