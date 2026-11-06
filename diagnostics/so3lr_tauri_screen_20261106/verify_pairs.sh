#!/bin/bash
set -euo pipefail
ROOT=$(git rev-parse --show-toplevel)
exec "$ROOT/workspaces/so3lr_tauri_screen_20261106/software/env/bin/python" "$ROOT/diagnostics/so3lr_tauri_screen_20261106/verify_pair_records.py"
