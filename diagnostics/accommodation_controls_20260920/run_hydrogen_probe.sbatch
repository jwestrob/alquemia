#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J H-minimize-probe
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -euo pipefail
PROBE_ROOT=${1:?absolute repository}
PROBE_RUN=${2:?absolute run directory}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
PROBE_SAMPLE=$((4-SLURM_ARRAY_TASK_ID))
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python "$PROBE_RUN/hydrogen_probe.py" --root "$PROBE_ROOT" --sample "$PROBE_SAMPLE" --output "$PROBE_RUN/sample_${PROBE_SAMPLE}"
