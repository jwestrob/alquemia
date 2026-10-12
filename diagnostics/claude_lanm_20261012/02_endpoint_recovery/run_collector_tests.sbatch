#!/bin/bash
#SBATCH -p standard
#SBATCH -J Nikasha-LanM-collector-tests
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -w node-64-768g-2
#SBATCH -o /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/claude_lanm_20261012/recovery_02_v1/collector_tests_%j.out

set -euo pipefail
if [[ $# -ne 1 ]]; then
    echo 'Usage: run_collector_tests.sbatch ABSOLUTE_FRESH_TEST_DIRECTORY' >&2
    exit 2
fi
TASK_GATE=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/claude_lanm_20261012/worktree_02/diagnostics/claude_lanm_20261012/02_endpoint_recovery/collector_test_gate.py
python3 "$TASK_GATE" run --directory "$1"
