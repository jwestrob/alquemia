#!/bin/bash
# Saved-array reproduction only; no new ML, DFT or MD calls.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
task_replay=$(mktemp -d "workspaces/so3lr_swap_spectator_20261106/replay_$(date -u +%Y%m%dT%H%M%SZ)_XXXXXX")
task_job=$(sbatch --parsable -N1 -n1 -c1 --qos=standard_expand --job-name=tauri_swap_saved_replay \
  --output="$task_replay/analysis_%j.out" --error="$task_replay/analysis_%j.err" \
  diagnostics/so3lr_swap_spectator_20261106/analyze.sh --output-dir "$PWD/$task_replay")
printf 'Saved-data replay job %s; outputs: %s\n' "$task_job" "$task_replay"
