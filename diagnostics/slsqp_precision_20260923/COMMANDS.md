# One frozen numerical precision pilot

Use the actual snapshots; do not submit a second attempt. GPU job1210561 owns
all eight searches and up to eight cross-MACE cells. Its wrapper then prepares
the finite native solvent manifest. The CPU wrapper executes that actual
manifest once and writes the final pool collection, retaining failures.

```bash
precision_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
precision_python=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
precision_run="$precision_root/workspaces/slsqp_precision_20260923"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$precision_python" "$precision_run/proposals_v1/implementation/slsqp_precision.py" validate \
  --manifest "$precision_run/proposals_v1/manifest.json"
"$precision_python" "$precision_run/pool_v1/implementation/slsqp_precision.py" validate_pool \
  --manifest "$precision_run/pool_v1/manifest.json"
"$precision_python" "$precision_run/pool_v1/implementation/slsqp_precision.py" compare \
  --collection "$precision_run/pool_v1/collection_final.json" \
  --output "$precision_run/COMPARISON_v1.json"
```

The comparison is write-once: read its completed file if it exists. No molecular
calls occur in validation/comparison. No reference is recalibrated. Independent
source preparation is fully parameterized by `prepare --help`; original pinned
manifests and actual submission receipts are in this experiment's workspace.

```bash
"$precision_python" -m unittest discover -s "$precision_root/tests" \
  -p test_slsqp_precision.py -v
```
