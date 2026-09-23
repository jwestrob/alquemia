#!/bin/bash
# One declared batch; preparation followed by the existing unchanged executor.
# Submit with 32 tasks x 1 CPU, one H200 and 200000 MiB as documented in the SOP.
set -euo pipefail
PLM_REQUEST=$1
case "$PLM_REQUEST" in /*) ;; *) exit 2;; esac
PLM_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
PLM_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
PLM_RUN=$(dirname "$PLM_REQUEST")/prepared
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
"$PLM_PY" "$PLM_ROOT/scripts/pqq_plm_prepare.py" prepare --request "$PLM_REQUEST" --output "$PLM_RUN"
if [ ! -f "$PLM_RUN/scoring/plan.json" ]; then
  echo "No complete supported triples. See $PLM_RUN/READY.json"
  exit 3
fi
bash "$PLM_ROOT/diagnostics/pqq_three_source_envelope_execution_20260923/run.sbatch" "$PLM_RUN/scoring/plan.json"
