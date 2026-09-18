#!/bin/bash
#SBATCH -p standard
#SBATCH -J alquemia-CHELPG-identity
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_chelpg_density_identity.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
IDENTITY_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
IDENTITY_ROOT=$(dirname "$1")
IDENTITY_DRIVER=$IDENTITY_ROOT/implementation/mace_chelpg_sampling.py
"$IDENTITY_PY" "$IDENTITY_DRIVER" dry-run-identity --manifest "$1"
"$IDENTITY_PY" "$IDENTITY_DRIVER" execute-identity --manifest "$1" --output "$IDENTITY_ROOT/collection_job_$SLURM_JOB_ID.json"
