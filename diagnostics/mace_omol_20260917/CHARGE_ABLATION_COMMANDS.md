# Charge-feature ablation operations

Protocol: `mace_omol_intact_charge_feature_ablation_descriptor_v1`.
These are modified model descriptors, not quantum endpoint energies.
Read CHARGE_ABLATION_PLAN.md and the active goal before changing the experiment.

## Existing run

Job 1200828, 42 forwards, one A5000 / 16 CPUs / 64474 MiB. Core numerical checks
precede whole-chain work; alpha numerical checks precede other structures.
No forces, reference recalibration, geometry search or production change.
The scheduler QOS is not a project compute budget. Do not duplicate a live job.

Use these explicit paths on biotite:

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
MACE_PY="$PWD/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python"
MACE_RUN="$PWD/workspaces/mace_omol_20260917/charge_ablation_development_v2"
"$MACE_PY" "$MACE_RUN/implementation/mace_hybrid.py" dry-run --manifest "$MACE_RUN/manifest.json"
squeue -j 1200828 -o '%.10i %.30j %.10T %.10M %.6C %.28R'
```

The exact submission argv and response are in `$MACE_RUN/submission.json`.
The existing wrapper collects actual results into `collection_job_1200828.json`.
A missing/failed endpoint remains unavailable. Rerunning the existing executor
under a fresh allocation uses verified completed tasks and preserves all attempts;
never delete its lock, overwrite its manifest, or run two executors for it.

## Collect, compare and report

After the job finishes, from the same directory and explicit variables above:

```bash
"$MACE_PY" "$MACE_RUN/implementation/mace_omol_ablation_run.py" report \
  --manifest "$MACE_RUN/manifest.json" \
  --output "$PWD/workspaces/mace_omol_20260917/charge_ablation_report_v1"
```

This validates the manifest and actual receipts, computes the declared paired
contrasts and distant-sodium consistency check, and writes unrounded values,
source endpoints, missing statuses, and Markdown. The output directory must be
new. Reports label modified outputs as kcal-equivalent model units. The baseline
and original native MACE values remain separate; no calibrated band exists.
Only `canonical_extension_permitted: true` permits the conditional 25-case run.

## Preparation and component audit

The component check is pinned under `charge_ablation_component_v3/`: 201 charge
categories, exact zero-charge/original-spin projection, unchanged parameters,
zero molecular energy calls. V1/v2 failed technical checks remain preserved.

Preparation uses `mace_omol_ablation_run.py prepare` with explicit
`--core-collection`, `--numerical-collection`, `--source-report`,
`--spectator-collection`, `--component`, `--agreement`, and `--output` arguments.
The exact source paths and hashes are in the prepared manifest. It contains 8
core checks, 14 alpha checks, 16 other primaries and 4 GGR sodium states.
All coordinates, charges, electron counts, waters, spin and source mappings
match their pinned prior tasks. The changed model feature has its own cache key.
V1 preparation was never executed; V2 corrects receipt validation for MACE's
atom-expanded charge embedding. The scientific model and settings are identical.

## Validation to date

Four real-fixture preflight tests pass (13.062 s); six existing native OMOL
regressions pass (14.587 s). No fabricated successful model output is used.
Scientific numerical and prediction checks belong to job 1200828 and must be
reported separately from these parser/state/cache tests. Do not infer physical
validity from feature invariance: it is imposed by the ablation architecture.


## Conditional canonical extension

Implementation is available as `scripts/mace_omol_ablation_panel.py`. It rejects
an incomplete, failed, changed or native-model development qualification, and
pins the strict v2 preparation (27 supported, 1KB0 unsupported). The planned
extension contains 100 new calibration forwards and eight exact descriptor
crystal reuses, with all 28 evidence rows retained.

Only after job 1200828 finishes and its full development gate passes:

```bash
"$MACE_PY" scripts/mace_omol_ablation_panel.py prepare \
  --preparation "$PWD/workspaces/mace_omol_20260917/intact_panel_prepared_v2/preparation_manifest.json" \
  --development-collection "$MACE_RUN/collection_job_1200828.json" \
  --agreement "$PWD/diagnostics/mace_omol_20260917/CHARGE_ABLATION_PLAN.md" \
  --output "$PWD/workspaces/mace_omol_20260917/charge_ablation_canonical_v1"
```

Preparation verifies the gate itself before creating a new manifest. Do not
supply a partial progress collection or change the declared criteria. The
canonical runner reuses the same executor and GPU allocation convention. A
new model only receives bands if all 25 calibration cases separate under the
frozen rule; the unavailable 1KB0 never becomes an implicit successful transfer.
Eight combined real-fixture/actual-forward guard tests pass (57.308 s), including
rejection of a native qualification and incomplete ablation collection.


Job **1200830 submitted** after the gate passed; both preparation and frozen
dry-run passed. Manifest SHA:
`2415b4dacd755530dec888017d76afe8895708213a36ef09497c7264e341476f`.
After its completed collection is present:

```bash
MACE_CAN="$PWD/workspaces/mace_omol_20260917/charge_ablation_canonical_v1"
"$MACE_PY" "$MACE_CAN/implementation/mace_omol_ablation_panel.py" report \
  --collection "$MACE_CAN/collection_job_1200830.json" \
  --output "$PWD/workspaces/mace_omol_20260917/charge_ablation_canonical_report_v1"
```

Use the frozen implementation so later code changes cannot silently change this
experiment. A previously created report directory is never overwritten.
