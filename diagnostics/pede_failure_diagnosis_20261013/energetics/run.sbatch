#!/bin/bash
#SBATCH --job-name=Nikasha-PedE-components
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
set -euo pipefail
: "${SLURM_JOB_ID:?Actual Slurm allocation required}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
config=${1:?Pass the exact frozen CONFIG.json}
python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pede_failure_diagnosis_20261013/energetics/ANALYZE.py --config "$config"
