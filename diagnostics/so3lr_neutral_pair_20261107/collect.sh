#!/bin/bash
set -euo pipefail
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TAURI_REPORT_DIR="${TAURI_REPORT_DIR:-$PWD/workspaces/so3lr_neutral_pair_20261107/saved_checks_${SLURM_JOB_ID}}"
mkdir -p "$TAURI_REPORT_DIR"
test ! -e "$TAURI_REPORT_DIR/CONTROL_SUMMARY.json"
cp diagnostics/so3lr_neutral_pair_20261107/{collect.py,audit_provenance.py,collect.sh,CHECKS.md,util.py,trajectory.py,OBSERVABLES.md} "$TAURI_REPORT_DIR/"
python diagnostics/so3lr_neutral_pair_20261107/audit_provenance.py
/groups/banfield/users/jwestrob/envs/lnsel/bin/python diagnostics/so3lr_neutral_pair_20261107/collect.py
if test -f workspaces/so3lr_neutral_pair_20261107/hydrated_preparation_v1/RESULT.json; then
 workspaces/so3lr_tauri_screen_20261106/software/env/bin/python diagnostics/so3lr_neutral_pair_20261107/trajectory.py --output-dir "$TAURI_REPORT_DIR/original_initial"
 workspaces/so3lr_tauri_screen_20261106/software/env/bin/python diagnostics/so3lr_neutral_pair_20261107/trajectory.py --preparation workspaces/so3lr_neutral_pair_20261107/hydrated_preparation_v1/RESULT.json --output-dir "$TAURI_REPORT_DIR/hydrated_initial"
 python - <<'PY'
import json,os
from pathlib import Path
w=Path('workspaces/so3lr_neutral_pair_20261107');out=Path(os.environ['TAURI_REPORT_DIR'])
for new,old in [('original_initial','initial_observations_53123'),('hydrated_initial','hydrated_observations_53157')]:
 a=json.loads((out/new/'RESULT.json').read_text());b=json.loads((w/old/'RESULT.json').read_text())
 assert a['rows']==b['rows'] and a['proton_transfer_events']==b['proton_transfer_events']
(out/'INITIAL_REPLAY_CHECKS.json').write_text(json.dumps(dict(passed=True,exact_real_records=24,molecular_calls=0),indent=2)+'\n')
PY
fi
