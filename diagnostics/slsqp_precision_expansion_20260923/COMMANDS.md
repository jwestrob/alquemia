# Precision34: finite operations and actual jobs

GPU1210813 owns all60 new searches and60 cross-MACE cells. CPU1211010 owns
exactly240 native GFN2 tasks in the generated pool manifest. Four completed
precision pools are reused. No q0, DFT or additional starts are authorized here.
Do not resubmit either recorded job.

```bash
precision_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
precision_python=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
precision_run="$precision_root/workspaces/slsqp_precision_expansion_20260923"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$precision_python" "$precision_run/proposals_v1/implementation/slsqp_precision_expansion.py" validate \
  --manifest "$precision_run/proposals_v1/manifest.json"
"$precision_python" "$precision_run/pool_v1/implementation/slsqp_precision_expansion.py" validate_pool \
  --manifest "$precision_run/pool_v1/manifest.json"
```

CPU wrapper collects actual outputs to pool_v1/collection_final.json, preserving
failure states. Report-only comparison and independent canonical25 calibration
(use once; outputs are immutable):

```bash
"$precision_python" "$precision_root/scripts/slsqp_precision_expansion_report.py" \
  --collection "$precision_run/pool_v1/collection_final.json" \
  --output "$precision_run/COMPARISON_v1.json" \
  --reference-output "$precision_run/REFERENCE_v1.json"
"$precision_python" -m unittest discover -s "$precision_root/tests" -p test_slsqp_precision_expansion.py -v
```

The reporting implementation is separate from the immutable molecular snapshot.
Original union225 settings/results remain unchanged. Both old-reference calls
and the separate new-reference calls must remain visible. Two old optimizer
failures have no successful-old pool/energy comparison; their final unsuccessful
iterates are diagnostic only. The archived repaired-H MMOL1770 source is a
separate experiment and deliberately retains the original optimizer precision.

## Completed outcome

Both jobs are terminal and both write-once report outputs exist. Read REPORT.md
and RESULT.json:34/34 calls retained and both old search failures recovered,
but strict numerical equivalence is27/32. Do not rerun completed chemistry or
overwrite the frozen reference. The read-only tests now check final actual
results and failed gates. Any full225 transfer is a separate prepared candidate,
with its rank policy awaiting independent qualification; no such jobs launched.
