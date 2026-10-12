#!/bin/bash
# Corrected replacement for diagnostics/metal_environment_response_20260926/run_lady_embedded.sbatch
# (FND-1, review_01_v1/FINDINGS.json): that wrapper hardcoded a 2026-09-28 shutdown
# cutoff that now unconditionally exits 75, and requested --ntasks=344 for a
# 112-rank x 2-concurrent-task layout this standard-QOS singleton-origin
# continuation does not use (QOS caps standard jobs at 12 CPUs; see
# workspaces/orca_parameter_audit_20261012/slurm_resumption_20261012/SCHEDULER_POLICY.json).
# This wrapper has no date cutoff and derives --ntasks from the manifest's own
# execution_resources (mpi_ranks * concurrent_tasks) instead of a fixed number,
# so it cannot silently drift out of sync with a differently-sized manifest again.
#SBATCH --job-name=Nikasha-LaDy-singleton-origin
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue

set -euo pipefail
if [ "$#" -ne 1 ]; then
    echo 'Usage: run_lady_embedded_singleton.sbatch ABSOLUTE_MANIFEST' >&2
    exit 2
fi
FIELD_MANIFEST="$1"
FIELD_WORK=$(dirname "$FIELD_MANIFEST")
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FIELD_SCRIPT="$FIELD_WORK/implementation/metal_environment_lady_embedded.py"
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export PATH=/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

REQUIRED_NTASKS=$("$FIELD_PY" -c "
import json
m = json.load(open('$FIELD_MANIFEST'))
r = m['execution_resources']
print(r['mpi_ranks'] * r['concurrent_tasks'])
")
if [ "${SLURM_NTASKS:-0}" -ne "$REQUIRED_NTASKS" ]; then
    echo "Allocation mismatch: SLURM_NTASKS=${SLURM_NTASKS:-0} but manifest requires $REQUIRED_NTASKS (mpi_ranks*concurrent_tasks)." >&2
    exit 76
fi

finish_field_reference() {
    "$FIELD_PY" "$FIELD_SCRIPT" collect --manifest "$FIELD_MANIFEST" \
        --output "$FIELD_WORK/COLLECTION_${SLURM_JOB_ID}.json" || true
    sacct -j "$SLURM_JOB_ID" --parsable2 --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,ReqMem,MaxRSS,ExitCode \
        > "$FIELD_WORK/accounting_in_job_${SLURM_JOB_ID}.txt" || true
}
trap finish_field_reference EXIT
"$FIELD_PY" "$FIELD_SCRIPT" execute --manifest "$FIELD_MANIFEST" \
    --output "$FIELD_WORK/EXECUTION_${SLURM_JOB_ID}.json" --force-gate "/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/metal_environment_response_20260926/dy_frozen_f_direction_v1/FINAL_COLLECTION.json"
