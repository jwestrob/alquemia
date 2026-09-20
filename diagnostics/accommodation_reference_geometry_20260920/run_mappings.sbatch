#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J carboxylate-mapping-check
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
CARBOXYLATE_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
CARBOXYLATE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$CARBOXYLATE_PY" "$CARBOXYLATE_ROOT/workspaces/accommodation_reference_geometry_20260920/implementation_v1/accommodation_reference_geometry.py" mappings --inventory "$CARBOXYLATE_ROOT/diagnostics/accommodation_reference_geometry_20260920/INVENTORY.json" --output "$CARBOXYLATE_ROOT/workspaces/accommodation_reference_geometry_20260920/mappings_v1" --workers 64
