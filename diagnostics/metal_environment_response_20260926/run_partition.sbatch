#!/bin/bash
#SBATCH --job-name=Nikasha-field-partition
#SBATCH --partition=high-memory,min_8t
#SBATCH --nodes=1
#SBATCH --ntasks=344
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --no-requeue

set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run_force_checks.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
FIELD_MANIFEST="$1"
FIELD_WORK=$(dirname "$FIELD_MANIFEST")
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FIELD_SCRIPT="$FIELD_WORK/implementation/metal_environment_partition.py"
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# Worker/rank layout comes from the new manifest; no fixed 4x16 override.
# Normal scheduler priority: no dependency on another session and no queue guard.
# Monday shutdown: do not begin a molecular stage after the conservative cutoff.
if [ "$(date +%s)" -ge "$(TZ=America/Los_Angeles date -d '2026-09-28 10:00:00' +%s)" ]; then
    echo 'Deferred: shutdown cutoff reached; recover explicitly after restart.' >&2
    exit 75
fi
finish_field_reference() {
    "$FIELD_PY" "$FIELD_SCRIPT" collect --manifest "$FIELD_MANIFEST" \
        --output "$FIELD_WORK/COLLECTION_${SLURM_JOB_ID}.json" || true
    sacct -j "$SLURM_JOB_ID" --parsable2 --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,ReqMem,MaxRSS,ExitCode \
        > "$FIELD_WORK/accounting_in_job_${SLURM_JOB_ID}.txt" || true
}
trap finish_field_reference EXIT
"$FIELD_PY" "$FIELD_SCRIPT" execute --manifest "$FIELD_MANIFEST" \
    --output "$FIELD_WORK/EXECUTION_${SLURM_JOB_ID}.json"
