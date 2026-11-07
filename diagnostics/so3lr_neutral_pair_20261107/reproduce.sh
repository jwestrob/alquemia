#!/bin/bash
# Reanalyse saved real arrays in a fresh directory; no ML, DFT or MD calls.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
task_replay=$(mktemp -d "workspaces/so3lr_neutral_pair_20261107/replay_$(date -u +%Y%m%dT%H%M%SZ)_XXXXXX")
export TAURI_REPORT_DIR="$PWD/$task_replay"
task_job=$(sbatch --parsable -N1 -n1 -c1 --qos=standard_expand --job-name=tauri_neutral_saved_replay \
 --output="$task_replay/analysis_%j.out" --error="$task_replay/analysis_%j.err" \
 diagnostics/so3lr_neutral_pair_20261107/collect.sh)
printf 'Saved-data replay job %s; outputs: %s\n' "$task_job" "$task_replay"
