#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J fold-native-DFT
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=600G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
FOLD_DFT_MASTER=${1:?absolute master manifest}
FOLD_DFT_SHARD=${2:?absolute shard manifest}
case "$FOLD_DFT_MASTER:$FOLD_DFT_SHARD" in /*:/*) ;; *) exit 2;; esac
FOLD_DFT_RUN=$(dirname "$FOLD_DFT_MASTER")
FOLD_DFT_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FOLD_DFT_ORCA=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PATH="$FOLD_DFT_ORCA:$PATH"
export LD_LIBRARY_PATH="$FOLD_DFT_ORCA${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
set +e
"$FOLD_DFT_PY" "$FOLD_DFT_RUN/implementation/accommodation_fold_dft.py" execute --manifest "$FOLD_DFT_SHARD"
FOLD_DFT_EXEC_RC=$?
"$FOLD_DFT_PY" "$FOLD_DFT_RUN/implementation/accommodation_fold_dft.py" collect --manifest "$FOLD_DFT_MASTER" --output "$FOLD_DFT_RUN/collection_job_${SLURM_JOB_ID}.json" --final-if-complete
FOLD_DFT_COLLECT_RC=$?
if (( FOLD_DFT_EXEC_RC != 0 )); then exit "$FOLD_DFT_EXEC_RC"; fi
exit "$FOLD_DFT_COLLECT_RC"
