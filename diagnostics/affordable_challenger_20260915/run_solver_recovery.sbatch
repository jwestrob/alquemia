#!/bin/bash
#SBATCH -p standard,memory,high-memory
#SBATCH -J alquemia-apbs-recover
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/affordable_challenger_20260915/solver_recovery_%j.out
#SBATCH -e /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/affordable_challenger_20260915/solver_recovery_%j.err
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/affordable_solver.py \
  --pilot-manifest workspaces/affordable_challenger_20260915/pilot/pilot_manifest.json \
  --output workspaces/affordable_challenger_20260915/solver_recovery \
  --apbs workspaces/affordable_challenger_20260915/software/apbs-3.4.1/APBS-3.4.1.Linux/bin/apbs \
  --prior-receipt workspaces/affordable_challenger_20260915/solver_terminal_receipt.json
