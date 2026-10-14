#!/bin/bash
#SBATCH --job-name=exaa-packing-prep
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 OPENMM_CPU_THREADS=1
exec /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/exaa_source_packing_20261014/packing/PREPARE.py \
 execute --manifest "$1"
