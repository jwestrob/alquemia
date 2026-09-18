#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-conductor-discrimination
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --array=0-6
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_conductor_discrimination.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=64 OMP_MAX_ACTIVE_LEVELS=1 OMP_PROC_BIND=spread OMP_PLACES=cores
CONDUCTOR_BASE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
CONDUCTOR_ROOT=$(dirname "$1")
CONDUCTOR_DRIVER=$CONDUCTOR_ROOT/implementation/mace_conductor_discrimination.py
CONDUCTOR_PY=$("$CONDUCTOR_BASE_PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["python"]["path"])' "$1")
CONDUCTOR_CASE=$("$CONDUCTOR_BASE_PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["groups"][int(sys.argv[2])]["case_id"])' "$1" "$SLURM_ARRAY_TASK_ID")
/usr/bin/time -v -o "$CONDUCTOR_ROOT/job_${SLURM_JOB_ID}.resources" \
    "$CONDUCTOR_PY" "$CONDUCTOR_DRIVER" execute --manifest "$1" --case "$CONDUCTOR_CASE"
