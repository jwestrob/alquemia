# Reproduce the Tinker capability inspection

```bash
ALQUEMIA_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ALQUEMIA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$ALQUEMIA_PY" "$ALQUEMIA_ROOT/scripts/mace_tinker_capability.py" dry-run \
  --manifest "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/tinker_capability_v1/manifest.json"
"$ALQUEMIA_PY" -m json.tool \
  "$ALQUEMIA_ROOT/diagnostics/mace_omol_20260917/TINKER_CAPABILITY_RESULT.json"
```

Executed export command, recorded for reproducibility. Existing output is
immutable: this command refuses to overwrite it.

```bash
"$ALQUEMIA_PY" "$ALQUEMIA_ROOT/scripts/mace_tinker_capability.py" prepare \
  --framework-result "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/amoeba_capability_v1/result.json" \
  --parameters "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/tinker_sources_v1/tinker_git/params/amoebabio18.prm" \
  --build-receipt "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/tinker_software_v1/receipt_job_1201011.json" \
  --states "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/full_boundary_GB_v1/states" \
  --plan "$ALQUEMIA_ROOT/diagnostics/mace_omol_20260917/TINKER_CAPABILITY_PLAN.md" \
  --output "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/tinker_capability_v1"
```

Executed parameter-only native operation; no energies/forces are available:

```bash
"$ALQUEMIA_PY" "$ALQUEMIA_ROOT/scripts/mace_tinker_capability.py" execute \
  --manifest "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/tinker_capability_v1/manifest.json"
```

Partial completed task records are reused only with identical configuration.
Complete result files are not overwritten. Inspect the existing result instead
of repeating the operation. Frozen implementation copies are in the workspace.

Tests, already passed without new native runs:

```bash
cd "$ALQUEMIA_ROOT"
"$ALQUEMIA_PY" -m unittest discover -s tests -p test_mace_tinker_capability.py -v
```

Build1201011 was submitted with `build_tinker.sbatch` and absolute manifest
`workspaces/mace_omol_20260917/tinker_software_v1/manifest_v2.json` under the
root above. It reuses partial objects from failed1201010 and never changes the
upstream source. The manifest, build driver, CMake compatibility include,
compiler/link commands and receipts are all preserved there. No installation
step or environment activation is required for the pinned parameter executable.
