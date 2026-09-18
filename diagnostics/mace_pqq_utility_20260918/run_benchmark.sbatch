#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
# Resources and method are explicit submission arguments; no project time cap.
set -euo pipefail
BENCH_DRIVER=$1
BENCH_MANIFEST=$2
BENCH_METHOD=$3
: "${SLURM_JOB_ID:?}" "${SLURM_CPUS_PER_TASK:?}"
BENCH_ORCA=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg
export ORCA_PATH="$BENCH_ORCA"
export PATH="$BENCH_ORCA:$PATH"
export LD_LIBRARY_PATH="$BENCH_ORCA${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export PYTHONUNBUFFERED=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
BENCH_WORK=$(dirname "$BENCH_MANIFEST")
scontrol show job "$SLURM_JOB_ID"
lscpu
if [[ "$BENCH_METHOD" == MACE ]]; then
    nvidia-smi --query-gpu=uuid,name,memory.total,driver_version --format=csv
fi
/usr/bin/time -v -o "$BENCH_WORK/${BENCH_METHOD}_job_${SLURM_JOB_ID}.resources.txt" \
    "$BENCH_DRIVER" "$BENCH_WORK/implementation/benchmark.py" execute \
    --manifest "$BENCH_MANIFEST" --method "$BENCH_METHOD"
