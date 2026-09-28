#!/bin/bash
#SBATCH --job-name=Nikasha-LaDy-summary
#SBATCH --partition=standard-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
set -euo pipefail
RESEARCH_ROOT="$1"
SUMMARY_DIRECTORY="$2"
FIELD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$RESEARCH_ROOT"
mkdir "$SUMMARY_DIRECTORY"
COMPARE_DIR="$RESEARCH_ROOT/diagnostics/metal_environment_response_20260926/repaired_exchange_20260928"
sha256sum "$COMPARE_DIR/COMPARE.py" "$COMPARE_DIR/COMPARISON_WITH_CLASSICAL_v2.json" "$COMPARE_DIR/COMPARISON_RECOVERY_v1.json" > "$SUMMARY_DIRECTORY/analysis_hashes.txt"
# These independent tables preserve the failed primary cells and separately named recovery.
summary_status=0
"$FIELD_PY" "$COMPARE_DIR/COMPARE.py" --config "$COMPARE_DIR/COMPARISON_WITH_CLASSICAL_v2.json" --output "$SUMMARY_DIRECTORY/primary.json" > "$SUMMARY_DIRECTORY/primary.stdout" 2> "$SUMMARY_DIRECTORY/primary.stderr" || summary_status=1
"$FIELD_PY" "$COMPARE_DIR/COMPARE.py" --config "$COMPARE_DIR/COMPARISON_RECOVERY_v1.json" --output "$SUMMARY_DIRECTORY/recovery.json" > "$SUMMARY_DIRECTORY/recovery.stdout" 2> "$SUMMARY_DIRECTORY/recovery.stderr" || summary_status=1
sacct -j 1220323,1220324,1220332,1220333,1220334,1220335,1220336,1220337,1220342,1220343 --parsable2 --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,ReqMem,MaxRSS,ExitCode > "$SUMMARY_DIRECTORY/endpoint_accounting.txt"
exit "$summary_status"
