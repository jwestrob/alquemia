#!/bin/bash
#SBATCH -p standard
#SBATCH -J Nikasha-Ca-GFN2
#SBATCH --ntasks=64
#SBATCH --nodes=1
#SBATCH -w node-64-768g-2
#SBATCH -o workspaces/orca_parameter_audit_20261012/GFN2_slurm_%j.out

set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo 'Usage: GFN2_STAGE.sbatch ABSOLUTE_SCOPE ABSOLUTE_TERMINAL_MACE_COLLECTION ABSOLUTE_OUTPUT_DIRECTORY' >&2
    exit 2
fi
GFN_SCOPE="$1"
GFN_COLLECTION="$2"
GFN_DIRECTORY="$3"
GFN_REPO=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
cd "$GFN_REPO"
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
gfn_accounting() {
    if [[ -d "$GFN_DIRECTORY" ]]; then
        sacct -j "$SLURM_JOB_ID" --parsable2 \
            --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,ReqMem,MaxRSS,ExitCode \
            > "$GFN_DIRECTORY/SLURM_ACCOUNTING_${SLURM_JOB_ID}.txt" || true
    fi
}
trap gfn_accounting EXIT
python3 diagnostics/orca_parameter_audit_20261012/slurm_resumption/GFN2_STAGE.py \
    --scope "$GFN_SCOPE" --collection "$GFN_COLLECTION" --directory "$GFN_DIRECTORY"
