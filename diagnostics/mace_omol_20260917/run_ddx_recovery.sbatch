#!/bin/bash
#SBATCH -p standard,memory
#SBATCH -J alquemia-ddPCM-recover
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --exclude=node-112-1500g-[1-4],node-224-3t-[1-2]
set -euo pipefail
if [[ $# -ne 2 || $1 != /* || $2 == */* ]]; then
    echo 'Usage: run_ddx_recovery.sbatch ABSOLUTE_MANIFEST GROUP_ID' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=64
DDX_BASE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
DDX_ROOT=$(dirname "$1")
DDX_DRIVER=$DDX_ROOT/implementation/mace_ddx_recovery.py
"$DDX_BASE_PY" "$DDX_DRIVER" dry-run --manifest "$1"
DDX_PY=$("$DDX_BASE_PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["python"]["path"])' "$1")
/usr/bin/time -v -o "$DDX_ROOT/job_${SLURM_JOB_ID}.resources" \
    "$DDX_PY" "$DDX_DRIVER" execute --manifest "$1" --group "$2" --output "$DDX_ROOT/group_${2}_job_$SLURM_JOB_ID.json"
