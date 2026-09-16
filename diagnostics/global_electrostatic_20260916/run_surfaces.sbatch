#!/bin/bash
#SBATCH -p standard-shared
#SBATCH -J alquemia-tabi
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
set -euo pipefail
if [ "$#" -lt 1 ]; then
    echo 'Usage: run_surfaces.sbatch ABSOLUTE_CAMPAIGN_MANIFEST [TASK_ID ...]' >&2
    exit 2
fi
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FIELD_MANIFEST=$1
shift
if [ "$#" -gt 0 ]; then
    "$FIELD_PY" scripts/global_electrostatic.py execute-surfaces --manifest "$FIELD_MANIFEST" --tasks "$@"
else
    "$FIELD_PY" scripts/global_electrostatic.py execute-surfaces --manifest "$FIELD_MANIFEST"
fi
