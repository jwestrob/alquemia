#!/bin/bash
#SBATCH --account=banfield
#SBATCH --qos=standard_expand
#SBATCH --partition=standard
#SBATCH --job-name=lbt-reference-direction-prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Actual allocation required}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE METAL_ENV_WORKERS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
unset OMP_PLACES OMP_PROC_BIND
export PATH="/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH"
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11 "$@"
