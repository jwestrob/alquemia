#!/bin/bash
#SBATCH -p standard
#SBATCH -J fold-proposal-solvent
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
set -uo pipefail
PROPOSAL_MANIFEST=$1
PROPOSAL_SHARD=$2
case "$PROPOSAL_MANIFEST" in /*) ;; *) exit 2;; esac
case "$PROPOSAL_SHARD" in 0|1|2|3) ;; *) exit 2;; esac
PROPOSAL_RUN=$(dirname "$PROPOSAL_MANIFEST")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
PROPOSAL_PYTHON=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$PROPOSAL_PYTHON" "$PROPOSAL_RUN/implementation/affordable_workflow.py" execute --manifest "$PROPOSAL_RUN/proposal_GFN2/shard_$PROPOSAL_SHARD/manifest.json"
SOLVENT_EXIT=$?
"$PROPOSAL_PYTHON" "$PROPOSAL_RUN/implementation/accommodation_proposals.py" collect --manifest "$PROPOSAL_MANIFEST" --output "$PROPOSAL_RUN/after_shard_${PROPOSAL_SHARD}_${SLURM_JOB_ID}.json"
COLLECTION_EXIT=$?
if [ "$SOLVENT_EXIT" -ne 0 ]; then exit "$SOLVENT_EXIT"; fi
exit "$COLLECTION_EXIT"
