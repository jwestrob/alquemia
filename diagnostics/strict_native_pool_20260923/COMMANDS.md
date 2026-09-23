# Exact operations

Run from the repository root. No production/default changes.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RUN="$PWD/workspaces/strict_native_pool_20260923/run_v1"
# Already prepared and submitted as1211270; do not resubmit completed work.
"$PY" "$RUN/fresh/implementation/strict_native_pool.py" validate --manifest "$RUN/fresh/manifest.json"
"$PY" "$RUN/cold_seed/implementation/strict_native_pool.py" validate --manifest "$RUN/cold_seed/manifest.json"
# Collection outputs are written automatically once each branch terminates.
# Analysis uses a separate frozen analysis implementation; write-once outputs.
"$PY" "$RUN/analysis_v1/strict_native_pool_compare.py" \
 --fresh "$RUN/fresh/COLLECTION.json" --cold-seed "$RUN/cold_seed/COLLECTION.json" \
 --old-reference "$PWD/workspaces/slsqp_precision_expansion_20260923/REFERENCE_v1.json" \
 --output "$RUN/COMPARISON.json" --reference-output "$RUN/REFERENCES.json"
"$PY" -m unittest discover -s tests -p test_strict_native_pool.py -v
```

Protocol `native_GFN2_precision32_TolE1e10_fresh_vs_cold_seed_v1`.
A uses fresh SAD/NoAutostart. B reuses original loose-cold GBW+xtbw before one strict call. One32CPU/64GiB allocation executes the384fresh calls followed by374newseeded calls; ten prior strict seeded results are explicitly reused. Reference/score fields remain separate; no other branch substitutes for failure.
