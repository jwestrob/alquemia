#!/bin/bash
#SBATCH --job-name=alquemia_khoury_masked
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=64474M
#SBATCH --time=7-00:00:00
set -euo pipefail
if [[ $# != 3 ]]; then
  echo 'Usage: sbatch run_khoury_benchmark.sbatch DRIVER_PYTHON MANIFEST_LIST PINNED_RUN_PILOT' >&2
  exit 2
fi
KHOURY_PY=$1
KHOURY_LIST=$2
KHOURY_RUNNER=$3
export MACE_MIN_MEMORY_MIB=64474
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
while IFS= read -r KHOURY_MANIFEST; do
  [[ -n "$KHOURY_MANIFEST" ]] || continue
  bash "$KHOURY_RUNNER" "$KHOURY_PY" "$KHOURY_MANIFEST" native || true
done < "$KHOURY_LIST"
