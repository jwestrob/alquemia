#!/bin/bash
# Rebuild saved analyses by default; --inference also repeats the fixed ML panel.
# Never reruns ORCA. Original numerical artifacts and reports remain intact.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
mode=${1:---saved}
[[ "$mode" == --saved || "$mode" == --inference ]] || { echo 'Use --saved or --inference' >&2; exit 2; }
task_dir=diagnostics/so3lr_tauri_screen_20261106
task_work=workspaces/so3lr_tauri_screen_20261106
task_replay=$(mktemp -d "$task_work/reproduce_$(date -u +%Y%m%dT%H%M%SZ)_XXXXXX")
export TAURI_RESULTS_DIR="$PWD/$task_replay/report"
mkdir -p "$TAURI_RESULTS_DIR"
task_result="$PWD/$task_work/panel_cpu_v1/RESULT.json"
task_dependency=''
if [[ "$mode" == --inference ]]; then
  task_attempt="$(basename "$task_replay")"
  task_ml=$(sbatch --parsable -N1 -n1 -c12 --qos=standard_expand --job-name=tauri_replay \
    --output="$task_replay/inference_%j.out" --error="$task_replay/inference_%j.err" \
    "$task_dir/cpu.sh" --attempt "$task_attempt")
  task_result="$PWD/$task_work/panel_cpu_${task_attempt}/RESULT.json"
  task_dependency=":$task_ml"
fi
task_ref=$(sbatch --parsable -N1 -n1 -c1 --qos=standard_expand --job-name=tauri_ref_replay \
  --output="$task_replay/reference_%j.out" --error="$task_replay/reference_%j.err" \
  "$task_dir/collect_references.sh")
task_pairs=$(sbatch --parsable -N1 -n1 -c1 --qos=standard_expand --job-name=tauri_pair_replay \
  --output="$task_replay/pair_%j.out" --error="$task_replay/pair_%j.err" "$task_dir/verify_pairs.sh")
task_analysis=$(sbatch --parsable -N1 -n1 -c1 --qos=standard_expand --job-name=tauri_analysis_replay \
  --dependency="afterok:$task_ref:$task_pairs$task_dependency" \
  --output="$task_replay/analysis_%j.out" --error="$task_replay/analysis_%j.err" \
  "$task_dir/analyze.sh" --result "$task_result")
printf 'Queued reference=%s pairs=%s analysis=%s; outputs: %s\n' "$task_ref" "$task_pairs" "$task_analysis" "$task_replay"
