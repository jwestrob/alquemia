#!/bin/bash
#SBATCH --job-name=nikasha-twoCa-prep
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=6
set -euo pipefail
: "${MANIFEST:?Set exact frozen manifest path}"
: "${SLURM_JOB_ID:?This preparation requires Slurm}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/ca_coverage_recovery_20261012/spectator_preparation/PREPARE.py \
  execute --manifest "$MANIFEST" --workers 5
