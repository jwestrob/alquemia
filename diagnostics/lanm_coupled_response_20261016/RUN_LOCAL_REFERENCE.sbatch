#!/bin/bash
#SBATCH --account=banfield
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --job-name=lbt48_reference
set -euo pipefail
if (( $# < 4 || $# > 5 )); then
    echo 'Usage: sbatch RUN_LOCAL_REFERENCE.sbatch ABS_ADAPTER CONFIG FRESH_OUTPUT scout|remaining|preflight [REMAINING_TASK]' >&2
    exit 2
fi
adapter=$1
config=$2
output=$3
stage=$4
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE METAL_ENV_WORKERS=1 PYTHONUNBUFFERED=1
unset OMP_PLACES OMP_PROC_BIND
export PATH="/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH"
args=(--stage "$stage" --config "$config" --output "$output")
if (( $# == 5 )); then args+=(--task "$5"); fi
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11 "$adapter" "${args[@]}"
