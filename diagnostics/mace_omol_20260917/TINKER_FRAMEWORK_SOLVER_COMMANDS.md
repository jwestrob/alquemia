# Native framework controls: runnable operations

All12energies completed as1201015; do not duplicate execution. Output creation
is exclusive, with matching completed task reuse on an explicitly needed restart.

```bash
ALQUEMIA_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ALQUEMIA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
TINKER_WORK=$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/tinker_framework_solver_v1
"$ALQUEMIA_PY" "$TINKER_WORK/implementation/mace_tinker_framework_solver.py" dry-run --manifest "$TINKER_WORK/manifest.json"
"$ALQUEMIA_PY" -m json.tool "$TINKER_WORK/collection_job_1201015.json"

# Recollection only: no native energy calculation; exclusive new output.
"$ALQUEMIA_PY" "$TINKER_WORK/implementation/mace_tinker_framework_solver.py" collect \
  --manifest "$TINKER_WORK/manifest.json" --output "$TINKER_WORK/recollection_v1.json"
cd "$ALQUEMIA_ROOT"
"$ALQUEMIA_PY" -m unittest discover -s tests -p test_mace_tinker_framework_solver.py -v
```

Full preparation replay, if needed, performs12native parameter initializations
but no energies. It refuses to overwrite existing output.

```bash
"$ALQUEMIA_PY" "$ALQUEMIA_ROOT/scripts/mace_tinker_framework_solver.py" prepare \
 --parent-manifest "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/tinker_capability_v1/manifest.json" \
 --software "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/tinker_framework_solver_software_v1/receipt.json" \
 --plan "$ALQUEMIA_ROOT/diagnostics/mace_omol_20260917/TINKER_FRAMEWORK_SOLVER_PLAN.md" \
 --output "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/tinker_framework_solver_replay_v1"
# Explicit execution of that new replay, only if needed:
# sbatch "$ALQUEMIA_ROOT/diagnostics/mace_omol_20260917/run_tinker_framework.sbatch" \
#  "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/tinker_framework_solver_replay_v1/manifest.json"
```
