#!/bin/bash
#SBATCH --job-name=alquemia_mace_pilot
#SBATCH --partition=gpu_h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --gres=gpu:1
#SBATCH --mem=257962M
#SBATCH --time=7-00:00:00
set -euo pipefail
if [[ $# -lt 3 ]]; then
  echo 'Usage: sbatch run_pilot.sbatch PYTHON MANIFEST MEMORY_MODE [TASK_ID ...]' >&2
  exit 2
fi
MACE_PY=$1
MACE_MANIFEST=$2
MACE_MEMORY_MODE=$3
shift 3
: "${SLURM_JOB_ID:?}" "${SLURM_CPUS_PER_TASK:?}"
MACE_SHARES=${MACE_GPU_SHARES:-1}
MACE_EXPECTED_MEMORY=$((2063701 * MACE_SHARES / 8))
if [[ "${SLURM_MEM_PER_NODE:-0}" -ne "$MACE_EXPECTED_MEMORY" ]]; then
  echo "Allocation mismatch: expected ${MACE_EXPECTED_MEMORY} MiB, received ${SLURM_MEM_PER_NODE:-missing}; no inference started." >&2
  exit 3
fi
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export PYTHONUNBUFFERED=1
MACE_RUNNER=$(dirname "$MACE_MANIFEST")/implementation/mace_hybrid.py
finish_collection() {
  MACE_EXIT=$?
  trap - EXIT
  "$MACE_PY" "$MACE_RUNNER" collect --manifest "$MACE_MANIFEST" \
    --output "$(dirname "$MACE_MANIFEST")/collection_job_${SLURM_JOB_ID}.json" || true
  exit "$MACE_EXIT"
}
trap finish_collection EXIT
MACE_TASK_ARGS=()
for MACE_TASK in "$@"; do MACE_TASK_ARGS+=(--task-id "$MACE_TASK"); done
nvidia-smi --query-gpu=index,uuid,name,memory.total,driver_version --format=csv
head -n 3 /proc/meminfo
scontrol show job "$SLURM_JOB_ID"
"$MACE_PY" "$MACE_RUNNER" execute --manifest "$MACE_MANIFEST" --memory-mode "$MACE_MEMORY_MODE" "${MACE_TASK_ARGS[@]}"
