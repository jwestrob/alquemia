#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-collective-scaffold
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
SCAFFOLD_MANIFEST=$1
case "$SCAFFOLD_MANIFEST" in /*) ;; *) exit 2;; esac
SCAFFOLD_RUN=$(dirname "$SCAFFOLD_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
SCAFFOLD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$SCAFFOLD_PY" "$SCAFFOLD_RUN/implementation/scaffold_accommodation.py" execute --manifest "$SCAFFOLD_MANIFEST"
SCAFFOLD_EXIT=$?
"$SCAFFOLD_PY" "$SCAFFOLD_RUN/implementation/scaffold_accommodation.py" collect --manifest "$SCAFFOLD_MANIFEST" --output "$SCAFFOLD_RUN/collection_${SLURM_JOB_ID}.json"
SCAFFOLD_COLLECT_EXIT=$?
if [ "$SCAFFOLD_EXIT" -ne 0 ]; then exit "$SCAFFOLD_EXIT"; fi
exit "$SCAFFOLD_COLLECT_EXIT"
