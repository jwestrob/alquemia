#!/bin/bash
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
set -euo pipefail
if [[ $# != 2 ]]; then
  echo 'Usage: sbatch run_cpu.sbatch PYTHON MANIFEST' >&2
  exit 2
fi
MACE_CPU_PY=$1
MACE_CPU_MANIFEST=$2
: "${SLURM_JOB_ID:?}" "${SLURM_CPUS_PER_TASK:?}"
if [[ "$SLURM_CPUS_PER_TASK" -ne 64 || "${SLURM_MEM_PER_NODE:-0}" -lt 131072 ]]; then
  echo 'Declared native CPU reference requires 64 CPUs and 128 GiB host RAM.' >&2
  exit 3
fi
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export PYTHONUNBUFFERED=1
MACE_CPU_RUNNER=$(dirname "$MACE_CPU_MANIFEST")/implementation/mace_hybrid.py
finish_cpu_collection() {
  MACE_CPU_EXIT=$?
  trap - EXIT
  "$MACE_CPU_PY" "$MACE_CPU_RUNNER" collect --manifest "$MACE_CPU_MANIFEST" \
    --output "$(dirname "$MACE_CPU_MANIFEST")/collection_job_${SLURM_JOB_ID}.json" || true
  exit "$MACE_CPU_EXIT"
}
trap finish_cpu_collection EXIT
lscpu
head -n 3 /proc/meminfo
scontrol show job "$SLURM_JOB_ID"
"$MACE_CPU_PY" "$MACE_CPU_RUNNER" execute --manifest "$MACE_CPU_MANIFEST" --memory-mode native
