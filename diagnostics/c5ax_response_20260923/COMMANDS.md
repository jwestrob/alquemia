# Read-only C5AX audit

From the repository root, choose a new write-once output path:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
scripts/audit_c5ax_response.py \
  --precision225 workspaces/slsqp_precision_transfer_20260923/COMPARISON_v1.json \
  --precision34 workspaces/slsqp_precision_expansion_20260923/COMPARISON_v1.json \
  --exclusions workspaces/preparation_coverage_20260923/EXCLUSION_AUDIT_v1.json \
  --strict32 workspaces/strict_native_pool_20260923/run_v1/COMPARISON.json \
  --strict225 workspaces/strict_native_transfer_20260923/run_v1/COMPARISON.json \
  --output workspaces/c5ax_response_20260923/AUDIT_replay.json
```

No executable scientific program is invoked. The source/result pins in the
output retain exact provenance; the replay requires the real archived artifacts.
