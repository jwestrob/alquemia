#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-reference-starts
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
STARTS_MANIFEST=$1
case "$STARTS_MANIFEST" in /*) ;; *) exit 2;; esac
STARTS_RUN=$(dirname "$STARTS_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
STARTS_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$STARTS_PY" "$STARTS_RUN/implementation/structure_informed_starts.py" execute --manifest "$STARTS_MANIFEST"
STARTS_EXIT=$?
"$STARTS_PY" "$STARTS_RUN/implementation/structure_informed_starts.py" collect --manifest "$STARTS_MANIFEST" --output "$STARTS_RUN/collection_${SLURM_JOB_ID}.json"
STARTS_COLLECT_EXIT=$?
if [ "$STARTS_EXIT" -ne 0 ]; then exit "$STARTS_EXIT"; fi
exit "$STARTS_COLLECT_EXIT"
