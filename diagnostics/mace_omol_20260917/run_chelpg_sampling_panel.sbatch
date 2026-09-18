#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-CHELPG-sampling
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_chelpg_sampling_panel.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
PANEL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
PANEL_ROOT=$(dirname "$1")
"$PANEL_PY" "$PANEL_ROOT/implementation/mace_chelpg_sampling_panel.py" dry-run --manifest "$1"
"$PANEL_PY" "$PANEL_ROOT/implementation/mace_chelpg_sampling_panel.py" execute --manifest "$1" --output "$PANEL_ROOT/population_job_$SLURM_JOB_ID.json"
