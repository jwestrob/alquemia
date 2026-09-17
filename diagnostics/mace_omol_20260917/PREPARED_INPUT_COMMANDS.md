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
