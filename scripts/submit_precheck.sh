#!/bin/bash
#SBATCH -p standard
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -J precheck_8dq2
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/logs/precheck_%j.out
#SBATCH -e /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/logs/precheck_%j.err

# No walltime per cluster discipline (handoff §"Hard constraints").
# --exclusive + $SLURM_CPUS_ON_NODE auto-fills whatever node we land on.

set -euo pipefail

ALCH_DIR=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
PYTHON=/home/jwestrob/miniconda3/envs/fep/bin/python3

cd "$ALCH_DIR"
mkdir -p logs

echo "host=$(hostname) cpus=$SLURM_CPUS_ON_NODE jobid=$SLURM_JOB_ID"
echo "started: $(date -Iseconds)"
echo "---"

"$PYTHON" -u scripts/precheck_8dq2.py \
    --workers "$SLURM_CPUS_ON_NODE" \
    "$@"

echo "---"
echo "finished: $(date -Iseconds)"
