#!/bin/bash
#SBATCH -p gpu
#SBATCH -J native-cell-recovery
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
RECOVERY_ROOT=$1
case "$RECOVERY_ROOT" in /*) ;; *) exit 2;; esac
RECOVERY_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
RECOVERY_STATUS=0
"$RECOVERY_PY" "$RECOVERY_ROOT/implementation/native_failed_cell_recovery.py" execute --manifest "$RECOVERY_ROOT/manifest.json" || RECOVERY_STATUS=1
"$RECOVERY_PY" "$RECOVERY_ROOT/implementation/native_failed_cell_recovery.py" collect --manifest "$RECOVERY_ROOT/manifest.json" --output "$RECOVERY_ROOT/COLLECTION.json" || RECOVERY_STATUS=1
exit "$RECOVERY_STATUS"
