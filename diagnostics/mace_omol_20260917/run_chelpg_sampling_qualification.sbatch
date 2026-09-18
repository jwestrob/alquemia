#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-CHELPG-replay
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 1 || $1 != /* ]]; then
    echo 'Usage: run_chelpg_sampling_qualification.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
REPLAY_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
REPLAY_ROOT=$(dirname "$1")
REPLAY_DRIVER=$REPLAY_ROOT/implementation/mace_chelpg_sampling.py
"$REPLAY_PY" "$REPLAY_DRIVER" dry-run --manifest "$1"
"$REPLAY_PY" "$REPLAY_DRIVER" execute --manifest "$1" --output "$REPLAY_ROOT/execution_$SLURM_JOB_ID.json"
