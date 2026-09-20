#!/bin/bash
#SBATCH -p gpu
#SBATCH -J mace-proposals
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200000M
#SBATCH --nodelist=node-224-2t-8gpu-1
set -uo pipefail
PROPOSAL_MANIFEST=$1
case "$PROPOSAL_MANIFEST" in /*) ;; *) exit 2;; esac
PROPOSAL_RUN=$(dirname "$PROPOSAL_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
PROPOSAL_PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$PROPOSAL_PYTHON" "$PROPOSAL_RUN/implementation/accommodation_proposals.py" propose --manifest "$PROPOSAL_MANIFEST"
PROPOSAL_EXIT=$?
"$PROPOSAL_PYTHON" "$PROPOSAL_RUN/implementation/accommodation_proposals.py" collect --manifest "$PROPOSAL_MANIFEST" --output "$PROPOSAL_RUN/after_proposals_${SLURM_JOB_ID}.json"
COLLECTION_EXIT=$?
if [ "$PROPOSAL_EXIT" -ne 0 ]; then exit "$PROPOSAL_EXIT"; fi
exit "$COLLECTION_EXIT"
