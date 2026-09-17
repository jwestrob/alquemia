# AMOEBA capability operations

No Context, energy, force, optimization or GPU execution is implemented here.
Framework preparation is separate from qualification of a full hybrid score.

Verify the frozen inputs, code and installed parameter sources:

```bash
ALQUEMIA_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ALQUEMIA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$ALQUEMIA_PY" "$ALQUEMIA_ROOT/scripts/mace_amoeba_capability.py" dry-run \
  --manifest "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/amoeba_capability_v1/manifest.json"
```

Report actual results (including any explicit preparation failures):

```bash
"$ALQUEMIA_PY" "$ALQUEMIA_ROOT/scripts/mace_amoeba_capability.py" report \
  --result "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/amoeba_capability_v1/result.json"
```

Executed preparation command, preserved for reproducibility; its output now
exists and rerunning it refuses to overwrite. Do not generate repeated new
preparations simply to reproduce the report.

```bash
"$ALQUEMIA_PY" "$ALQUEMIA_ROOT/scripts/mace_amoeba_capability.py" prepare \
  --inputs \
  "$ALQUEMIA_ROOT/workspaces/mace_global_benchmark_20260916/prepared_v1/GGR_1GLG/preparation.json" \
  "$ALQUEMIA_ROOT/workspaces/mace_global_benchmark_20260916/prepared_v1/ALPHA_1F6S/preparation.json" \
  "$ALQUEMIA_ROOT/workspaces/mace_global_benchmark_20260916/prepared_v1/ALPHA_6IP9/preparation.json" \
  --plan "$ALQUEMIA_ROOT/diagnostics/mace_omol_20260917/AMOEBA_CAPABILITY_PLAN.md" \
  --output "$ALQUEMIA_ROOT/workspaces/mace_omol_20260917/amoeba_capability_v1"
```

The preserved script and common helper under `amoeba_capability_v1/implementation/`
can reproduce that implementation against pinned inputs in a separately
declared output directory. CLI paths are explicit; preparation creates exactly
the three declared framework systems, never endpoint energies.

Tests, already passed:

```bash
cd "$ALQUEMIA_ROOT"
"$ALQUEMIA_PY" -m unittest discover -s tests -p test_mace_amoeba_capability.py -v
```

Missing real fixtures produce explicit skips. Corruption tests mutate copies
of actual preparations and never serve as scientific evidence.
