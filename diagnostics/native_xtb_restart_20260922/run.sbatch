#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-native-restart
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
RESTART_MANIFEST=$1
case "$RESTART_MANIFEST" in /*) ;; *) exit 2;; esac
RESTART_RUN=$(dirname "$RESTART_MANIFEST")
RESTART_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$RESTART_PY" "$RESTART_RUN/implementation/native_xtb_restart.py" execute --manifest "$RESTART_MANIFEST"
RESTART_EXIT=$?
"$RESTART_PY" "$RESTART_RUN/implementation/native_xtb_restart.py" collect --manifest "$RESTART_MANIFEST" --output "$RESTART_RUN/collection_${SLURM_JOB_ID}.json"
RESTART_COLLECT_EXIT=$?
if [ "$RESTART_EXIT" -ne 0 ]; then exit "$RESTART_EXIT"; fi
exit "$RESTART_COLLECT_EXIT"
