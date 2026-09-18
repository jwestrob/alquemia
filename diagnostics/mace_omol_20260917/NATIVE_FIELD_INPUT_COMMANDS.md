# Native supplied-field replay operations

Six solves completed as1201061; all nine checks pass. No rerun is needed.
Recollect into a new output path and run actual-artifact tests:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_input_v1/implementation/mace_native_field_input.py dry-run --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_input_v1/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_input_v1/implementation/mace_native_field_input.py collect --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_input_v1/manifest.json --output /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_input_v1/collection_replay_01.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/tests -p test_mace_native_field_input.py -v
```

The existing executor verifies successful receipts and reuses them. For a
partially completed manifest, the same explicit bounded submit operation is:

```bash
sbatch --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_input_v1/slurm_%j.out --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_input_v1/slurm_%j.err /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/mace_omol_20260917/run_native_field_input.sbatch /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/native_field_input_v1/manifest.json
```

Actual build commands/failed and successful sources are retained in
`native_field_input_software_v1` and `native_field_input_software_v2` beside
the workspace. `prepare --help` lists explicit inputs for a new output directory.
This replay accepts the four actual native fields only; do not relabel it a
quantum-source calculation or use its total as a calibrated metal score.
