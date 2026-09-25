#!/bin/bash
#SBATCH -p memory
#SBATCH -J nikasha-lanm-vacuum-seeded
#SBATCH --nodes=1
#SBATCH --ntasks=112
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jacobwestroberts@gmail.com
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
lanm_py=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
lanm_code=diagnostics/lanm_global_occupancy_20260923
lanm_out=workspaces/lanm_global_occupancy_20260923/seeded_vacuum_v1
read -r lanm_ranks lanm_memory lanm_maxcore < <("$lanm_py" "$lanm_code/native_allocation.py" --workers 2)
"$lanm_py" "$lanm_code/seed_vacuum_recovery.py" prepare --output "$lanm_out" --ranks "$lanm_ranks" --memory-mib "$lanm_memory" --maxcore-mb "$lanm_maxcore"
"$lanm_py" "$lanm_code/seed_vacuum_recovery.py" execute --manifest "$lanm_out/manifest.json"
