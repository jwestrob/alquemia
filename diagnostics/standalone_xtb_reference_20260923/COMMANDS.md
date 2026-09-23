# Canonical standalone reference operations

Job1210476 owns the300 new finite cells;36 actual pilot cells are reused.
Do not submit this manifest a second time. The wrapper writes collection_v1.json
when terminal. No225-fold calculation is included.

```bash
ref_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ref_python=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
ref_run="$ref_root/workspaces/standalone_xtb_reference_20260923/run_v1"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$ref_python" "$ref_run/implementation/standalone_xtb_reference.py" validate --manifest "$ref_run/manifest.json"
"$ref_python" "$ref_run/implementation/standalone_xtb_reference.py" calibrate \
  --collection "$ref_run/collection_v1.json" --output "$ref_run/REFERENCE_v1.json"
```

Calibration is write-once. If REFERENCE_v1.json already exists, read that record;
do not rerun or overwrite. It retains static/minimal/mathematical variants,
all28 rows, exact native scores,25canonical calibration records and3consumed
crystal-transfer decisions. It may legitimately record unsupported separation.

Read-only tests use actual prepared/computed artifacts:

```bash
"$ref_python" -m unittest discover -s "$ref_root/tests" \
  -p test_standalone_xtb_reference.py -v
```

Future source inventory is TRANSFER_INPUT_INVENTORY_v1.json one directory above
ref_run. It verifies208static sources/416actual MACE origins and206minimal pools,
retains the17preparation exclusions and original adaptive failures, and is
explicitly not an execution manifest or permission to launch225-fold scoring.
