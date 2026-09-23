#!/bin/bash
#SBATCH -p gpu
#SBATCH -J nikasha-cpcm
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
set -uo pipefail
CPCM_MANIFEST=$1
case "$CPCM_MANIFEST" in /*) ;; *) exit 2;; esac
CPCM_RUN=$(dirname "$CPCM_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
CPCM_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$CPCM_PY" "$CPCM_RUN/implementation/solvent_cpcm_pilot.py" execute --manifest "$CPCM_MANIFEST" >"$CPCM_RUN/execution_${SLURM_JOB_ID}.json"
CPCM_EXIT=$?
"$CPCM_PY" "$CPCM_RUN/implementation/solvent_cpcm_pilot.py" collect --manifest "$CPCM_MANIFEST" --output "$CPCM_RUN/collection_${SLURM_JOB_ID}.json"
CPCM_COLLECT_EXIT=$?
if [ "$CPCM_EXIT" -ne 0 ]; then exit "$CPCM_EXIT"; fi
exit "$CPCM_COLLECT_EXIT"
