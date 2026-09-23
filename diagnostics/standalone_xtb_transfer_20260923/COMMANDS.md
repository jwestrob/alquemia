# Frozen static225 standalone operations

Job1210508 owns the828 new cells. Four exact prior pilot cells are reused.
Do not resubmit the manifest. The wrapper writes the actual collection on
termination. No fresh adaptive, MACE, DFT or geometry search belongs to this run.

```bash
static_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
static_python=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
static_run="$static_root/workspaces/standalone_xtb_transfer_20260923/run_v1"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$static_python" "$static_run/implementation/standalone_xtb_transfer.py" validate --manifest "$static_run/manifest.json"
"$static_python" "$static_run/implementation/standalone_xtb_transfer.py" compare \
  --collection "$static_run/collection_v1.json" --output "$static_run/COMPARISON_v1.json"
```

Comparison is write-once: read its existing file after completion. It retains
all225 cases, strict groups, old native/DFT/adaptive methods, source preparation
failures and any fresh executable failures. Bands are the canonical reference's
unchanged static variant; no calibration occurs in this operation.

```bash
"$static_python" -m unittest discover -s "$static_root/tests" \
  -p test_standalone_xtb_transfer.py -v
```
