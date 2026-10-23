#!/bin/bash
#SBATCH --account=banfield
#SBATCH --qos=standard_expand
#SBATCH --partition=gpu_h200
#SBATCH --job-name=lady-tinker9-parity
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual allocation required}"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export LC_ALL=C
task_repo=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
exec /usr/bin/python3 "$task_repo/diagnostics/lanm_selectivity_free_energy_20261021/gpu_engine_qualification_v1/RUN.py" \
  --config "$task_repo/diagnostics/lanm_selectivity_free_energy_20261021/gpu_engine_qualification_v1/CONFIG_BOUND.json" \
  --output "$task_repo/workspaces/lanm_selectivity_free_energy_20261021/gpu_engine_qualification_v1/execution_v1"
