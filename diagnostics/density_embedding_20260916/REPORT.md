# Density and permanent-field diagnostic — running

Jacob's [autonomous improvement authorization and experiment scope](AGREEMENT.md)
is recorded before execution. Baseline/default and the failed global v1 record
remain unchanged. New protocol:
`native_r2scan3c_permanent_field_density_diagnostic_v1`.

Four actual 1H4I qm33/qm36 La/Ca states were prepared with exact source XYZ
files and byte-identical permanent fields within each metal pair. Three real-
artifact tests pass; runner dry-run and both batch-script syntax checks pass.
The code also preserves the existing field exclusions and native basis/ECP.

- Job **1199979**: four ORCA-vpot evaluations at all actual environmental charge
  positions, using copied saved vacuum wavefunctions/densities; four CPUs.
- Job **1199980**: four native vacuum r2SCAN-3c/MBIS calculations in those same
  permanent fields; four 16-rank workers on 64 CPUs.

No density-coupling or embedded result is reported yet. The [energy accounting](ACCOUNTING.md)
distinguishes exact-density coupling, monopole approximation, and electronic
response. Neither operation produces a calibrated affinity score.

Workspace: `workspaces/density_embedding_20260916/`. Explicit manifests:
`prepared_v1/{manifest,potential_manifest}.json`. Execution commands and actual
Slurm job IDs are retained in `{embedded,potentials}_submission.json`.
Implementation snapshots: `prepared_v1/implementation/`.

## Read-only collection

Run from the repository root after the corresponding jobs finish:

```bash
DENSITY_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
DENSITY_WORK=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/density_embedding_20260916
"$DENSITY_PY" scripts/density_embedding.py collect-potentials \
  --manifest "$DENSITY_WORK/prepared_v1/potential_manifest.json" \
  --output "$DENSITY_WORK/potentials_operator_v1.json"
"$DENSITY_PY" scripts/density_embedding.py collect-embedded \
  --manifest "$DENSITY_WORK/prepared_v1/manifest.json" \
  --potentials "$DENSITY_WORK/potentials_operator_v1.json" \
  --output "$DENSITY_WORK/embedded_operator_v1.json"
```

Writers refuse existing outputs. Missing/failed utility or quantum outputs
remain unavailable. This report will be updated from actual collected results.
