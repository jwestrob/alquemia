# Replay the completed native accounting report

Job1201055 already executed all18native operations. Six initial coordinate
validation errors were recovered from unchanged outputs, without reruns.
Use the pinned recovery below; the original job collection remains an immutable
record of the validation failure. Do not resubmit its old executor to recover
those already computed outputs.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_accounting_v1/recovery_implementation_v1/mace_native_field_accounting.py dry-run --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_accounting_v1/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_accounting_v1/recovery_implementation_v1/mace_native_field_accounting.py collect --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_accounting_v1/manifest.json --recovery /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_accounting_v1/recovery_pinned_v1/recovery.json --output /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_accounting_v1/collection_replay_01.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/tests -p test_mace_native_field_accounting.py -v
```

The output is new and must not already exist. Exact historical submit/build
commands and input hashes are in submission.json and the pinned software receipt.
For a new preparation, use the current driver `prepare` with explicit parent
manifest, software receipt, plan and a new output directory. The CLI `--help`
lists each argument. Both failed and recovered results are retained.
