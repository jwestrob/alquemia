#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-GGR-density-GK
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64000M
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_trial_gk_expansion_native.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=64
PANEL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
PANEL_ROOT=$(dirname "$1")
PANEL_DRIVER=$PANEL_ROOT/implementation/mace_density_panel.py
"$PANEL_PY" "$PANEL_DRIVER" dry-run --manifest "$1"
"$PANEL_PY" "$PANEL_DRIVER" execute --manifest "$1"
"$PANEL_PY" "$PANEL_DRIVER" collect --manifest "$1" --output "$PANEL_ROOT/collection_job_$SLURM_JOB_ID.json"
