#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TAURI_REPORT_DIR="${TAURI_REPORT_DIR:-$PWD/workspaces/so3lr_neutral_pair_20261107/saved_checks_${SLURM_JOB_ID}}"
mkdir -p "$TAURI_REPORT_DIR"
test ! -e "$TAURI_REPORT_DIR/CONTROL_SUMMARY.json"
cp diagnostics/so3lr_neutral_pair_20261107/{collect.py,audit_provenance.py,collect.sh,CHECKS.md,util.py} "$TAURI_REPORT_DIR/"
python diagnostics/so3lr_neutral_pair_20261107/audit_provenance.py
/groups/banfield/users/jwestrob/envs/lnsel/bin/python diagnostics/so3lr_neutral_pair_20261107/collect.py
