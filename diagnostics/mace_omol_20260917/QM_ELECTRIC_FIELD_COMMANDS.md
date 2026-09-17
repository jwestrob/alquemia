# Saved-density electric-field operations

Executed preparation is pinned under `qm_electric_field_v1`; execution1201017.
Check the live queue and its receipts before resuming; do not duplicate it.
These operations never launch DFT, fit charges or modify a baseline score.

```bash
ALQUEMIA_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ALQUEMIA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
FIELD_WORK=$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/qm_electric_field_v1
"$ALQUEMIA_PY" "$FIELD_WORK/implementation/mace_qm_field.py" dry-run --manifest "$FIELD_WORK/manifest.json"
"$ALQUEMIA_PY" "$FIELD_WORK/implementation/mace_qm_field.py" collect \
 --manifest "$FIELD_WORK/manifest.json" --output "$FIELD_WORK/recollection_v1"
cd "$ALQUEMIA_ROOT"
"$ALQUEMIA_PY" -m unittest discover -s tests -p test_mace_qm_field.py -v
```

Full input replay, no scientific calculation. Output creation is exclusive.

```bash
"$ALQUEMIA_PY" "$ALQUEMIA_ROOT/scripts/mace_qm_field.py" prepare \
 --charges "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/normalized_charge_report_v1/result.json" \
 --frameworks "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/amoeba_capability_v1/result.json" \
 --states "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/full_boundary_GB_v1/states" \
 --plan "$ALQUEMIA_ROOT/diagnostics/mace_omol_20260917/QM_ELECTRIC_FIELD_PLAN.md" \
 --output "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/qm_electric_field_replay_v1"
# Explicit new replay execution, only if needed:
# sbatch "$ALQUEMIA_ROOT/diagnostics/mace_omol_20260917/run_qm_field.sbatch" \
#  "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/qm_electric_field_replay_v1/manifest.json"
```

Execution reuses only an accepted identical task/manifest/executable receipt.
An interrupted/failed attempt remains visible and requires explicit
`execute --retry-failed` inside the declared8-worker allocation. This option
does not re-execute successful endpoints. Reports retain unavailable fields;
neither missing potentials nor a failed representation can become a zero
environmental correction. The `U0` diagnostic is never added to a score.
