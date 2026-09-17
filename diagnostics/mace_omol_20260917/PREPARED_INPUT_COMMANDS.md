# Prepared-input MACE descriptor operations

Opt-in research interface. Read PREPARED_INPUT_INTERFACE_REPORT.md and the active
CURRENT.md. No calibration/default promotion is implied.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
MACE_DRIVER=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
MACE_CASE="$PWD/workspaces/mace_omol_20260917/prepared_interface_ggr_fresh_v1"
"$MACE_DRIVER" "$MACE_CASE/implementation/mace_hybrid.py" dry-run \
  --manifest "$MACE_CASE/manifest.json"
```

That existing manifest has four unexecuted tasks. It is an interface example;
GGR already has exact completed descriptor results, so no duplicate submission
was made. The ordinary finite executor remains available for genuinely new
inputs. Use the recorded driver Python as the first argument to the existing
run_pilot.sbatch wrapper; its workers automatically use the separate pinned
MACE Python from the manifest's software record. Keep the existing one-A5000,
16CPU,64474MiB allocation and expandable-segments setting for this protocol.

The already completed no-inference PQQ interface report is:
`workspaces/mace_omol_20260917/prepared_interface_pqq_report_v1/result.json`.
It verifies/reuses four actual states and reports the exact original descriptor.
Classification is explicitly unavailable until compatible calibration is added.

For a new prepared input, the CLI requires explicit --preparation,
--development-collection, --agreement and --output paths. Optional
--reuse-collection requires actual matching receipts; it never accepts a baseline
cache. `audit --preparation` checks source replay without a model call.
The exact executable commands used for the current ready GGR example were:

```bash
"$MACE_DRIVER" scripts/mace_omol_prepared.py prepare \
  --preparation "$PWD/workspaces/mace_global_benchmark_20260916/prepared_v1/GGR_1GLG/preparation.json" \
  --development-collection "$PWD/workspaces/mace_omol_20260917/charge_ablation_development_v2/collection_job_1200828.json" \
  --agreement "$PWD/diagnostics/mace_omol_20260917/PREPARED_INPUT_INTERFACE_PLAN.md" \
  --output "$MACE_CASE"
```

That preparation already exists; this historical command intentionally refuses
to overwrite it. Select a new output directory when preparing genuinely new
inputs. No source edits are necessary to choose an input. After real execution,
the report operation takes explicit `--manifest` and `--output` paths and a
new output directory; unavailable endpoints remain unavailable.


## Optional exact two-call evaluation

The qualified factorization reference is:
`workspaces/mace_omol_20260917/factorization_report_v2/result.json`.
Supply it with `--factorization` during preparation to create only the two bound
states. The interface verifies actual saved readouts and the exact model identity;
it does not use a fitted or arbitrary offset. This reference is not a biological
calibration, aquo reference or physical ion energy. Four-state evaluation remains
available when the argument is omitted.

An executed interface check already exists at
`prepared_interface_pqq_two_call_v2/`, with two actual bound endpoint reuses and
zero new inference. Its complete report is
`prepared_interface_pqq_two_call_report_v2/result.json`. The optional path passed
three real-artifact tests; full-panel equivalence is pending the unchanged
canonical job. Any classification must use its own compatible numeric evaluation
record, as declared in FACTORIZATION_NUMERICAL_BANDS_ADDENDUM.md.


## Optional matching research calibration (2026-09-17)

The completed canonical test and189factorization checks now support an explicit
PQQ-only calibration. Read MASKED_CALIBRATION_REPORT.md. Production remains
unchanged;1KB0 remains unscorable and no generic affinity band is implied.
The following reports an already computed1H4I pair; it launches no model call:

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
MACE_DRIVER=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
MACE_WORK="$PWD/workspaces/mace_omol_20260917"
"$MACE_DRIVER" "$MACE_WORK/masked_calibration_source_v1/implementation/mace_omol_prepared.py" report \
  --manifest "$MACE_WORK/prepared_interface_pqq_two_call_v2/manifest.json" \
  --calibration "$MACE_WORK/masked_calibration_v1/reference.json" \
  --output "$MACE_WORK/prepared_interface_pqq_calibrated_review_v2"
```

The output path must be new; completed reports are never overwritten. The
calibration argument is optional. Without it, the interface returns only its
raw descriptor. Supplying this calibration for GGR returns the raw score and
an explicit out-of-scope decision, not an inherited PQQ label. Source input
must be an auditable whole-chain preparation, not arbitrary unprotonated XYZ.
