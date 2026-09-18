#!/bin/bash
#SBATCH -p standard
#SBATCH -J alquemia-CHELPG-density
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_chelpg_sampling_panel_identity.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PANEL_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
PANEL_ROOT=$(dirname "$1")
"$PANEL_PY" "$PANEL_ROOT/implementation/mace_chelpg_sampling_panel.py" dry-run-identity --manifest "$1"
"$PANEL_PY" "$PANEL_ROOT/implementation/mace_chelpg_sampling_panel.py" execute-identity --manifest "$1" --output "$PANEL_ROOT/collection_job_$SLURM_JOB_ID.json"
