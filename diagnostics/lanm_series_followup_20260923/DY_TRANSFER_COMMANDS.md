# Dy-source structural transfer — completed, no new calls

All three jobs1210367/1210368/1210372 are complete. Do not resubmit their molecular
tasks. Exact submission argument arrays are preserved in
DY_TRANSFER_SUBMISSIONS.json and DY_TRANSFER_CONTINUATION_SUBMISSION.json.

Read-only preflight:

```bash
lanm_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
lanm_python=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
lanm_transfer="$lanm_root/workspaces/lanm_series_followup_20260923/dy_transfer_v1"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$lanm_python" "$lanm_transfer/implementation/lanm_structural_transfer.py" validate \
  --manifest "$lanm_transfer/manifest.json"
```

The already executed collector command was:

```bash
"$lanm_python" "$lanm_transfer/implementation/lanm_structural_transfer.py" collect \
  --manifest "$lanm_transfer/manifest.json" \
  --mace-collection "$lanm_transfer/mace_results/result.json" \
  --continuation-manifest "$lanm_transfer/native_continuation/manifest.json" \
  --output "$lanm_transfer/final_collection.json"
```

Its output is immutable and already exists; read it instead of rerunning that
write-once operation. The new adapter retains six fresh Hans endpoints and
explicitly references the six previously completed Mex endpoints. There is no
missing-value baseline fallback or implicit source minimum.

Real-artifact mapping tests:

```bash
"$lanm_python" -m unittest discover -s "$lanm_root/tests" \
  -p test_lanm_structural_transfer.py -v
```

The original12endpoint pilot and its own seven tests remain separately preserved.
