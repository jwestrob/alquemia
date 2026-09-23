# Standalone xTB fixed-cell pilot

Job1210445 executes the immutable224-call manifest. Do not submit it again.
No new MACE/DFT/optimization is included. Exact argument arrays are in
SUBMISSION.json; the prepared runner snapshot remains unchanged.

```bash
xtb_root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
xtb_python=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
xtb_run="$xtb_root/workspaces/standalone_xtb_20260923/run_v1"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$xtb_python" "$xtb_run/implementation/standalone_xtb.py" validate --manifest "$xtb_run/manifest.json"
```

The wrapper automatically collects into collection_v1.json after execution.
Collection/parser fixes consume existing outputs only and use a new immutable
collection filename if needed. The corrected collection_v2.json is final; collection_v1 retains its parser failures. The following comparison has already completed:

```bash
"$xtb_python" "$xtb_root/scripts/standalone_xtb.py" compare \
  --collection "$xtb_run/collection_v2.json" --output "$xtb_run/COMPARISON_v1.json"
```

Outputs are write-once; read existing completed files instead of running a
command twice. Scientific execution is bounded by the exact manifest; no
accuracy escalation, retries, external/native mixed solvation or threshold fit.
