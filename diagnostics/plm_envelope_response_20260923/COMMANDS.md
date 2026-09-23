# Reproduce this read-only report

Run from the repository root. Choose unused output paths; products are write-once.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 scripts/plm_envelope_response.py analyze \
 --result workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RESULT_1211626.json \
 --output workspaces/plm_envelope_response_20260923/geometry_replay

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 scripts/plm_envelope_response.py plot \
 --geometry workspaces/plm_envelope_response_20260923/geometry_v1/GEOMETRY.json \
 --output workspaces/plm_envelope_response_20260923/figures_replay

mkdir -p workspaces/plm_envelope_response_20260923/tables_replay
python diagnostics/plm_envelope_response_20260923/export_tables.py \
 workspaces/plm_envelope_response_20260923/geometry_v1/GEOMETRY.json \
 workspaces/plm_envelope_response_20260923/tables_replay
```

No molecular evaluation, new optimization, score calibration or source selection occurs.
