#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-native-field-accounting
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64000M
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_native_field_accounting.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=64
FIELD_RUN_ROOT=$(dirname "$1")
FIELD_RUN_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FIELD_RUN_DRIVER=$FIELD_RUN_ROOT/implementation/mace_native_field_accounting.py
"$FIELD_RUN_PY" "$FIELD_RUN_DRIVER" dry-run --manifest "$1"
"$FIELD_RUN_PY" "$FIELD_RUN_DRIVER" execute --manifest "$1"
"$FIELD_RUN_PY" "$FIELD_RUN_DRIVER" collect --manifest "$1" \
    --output "$FIELD_RUN_ROOT/collection_job_${SLURM_JOB_ID}.json"
