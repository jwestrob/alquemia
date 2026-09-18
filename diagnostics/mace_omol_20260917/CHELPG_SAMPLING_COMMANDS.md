# Reproduce the completed charge-sampling comparison

These commands inspect/collect actual results without starting calculations.
Output directories must be new; existing scientific artifacts are never replaced.

```bash
ALQ_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ALQ_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
ALQ_WORK=$ALQ_ROOT/workspaces/mace_omol_20260917
export OPENBLAS_NUM_THREADS=1

"$ALQ_PY" "$ALQ_ROOT/scripts/mace_chelpg_sampling_panel.py" dry-run \
  --manifest "$ALQ_WORK/chelpg_sampling_panel_v1/manifest.json"

"$ALQ_PY" "$ALQ_ROOT/scripts/mace_chelpg_sampling_panel.py" collect \
  --manifest "$ALQ_WORK/chelpg_sampling_panel_v1/manifest.json" \
  --output "$ALQ_WORK/chelpg_sampling_review_collection_v1"

"$ALQ_PY" "$ALQ_WORK/chelpg_sampling_reporting_source_v1/mace_chelpg_sampling_report.py" \
  --population "$ALQ_WORK/chelpg_sampling_panel_v1/population_job_1201177.json" \
  --identity "$ALQ_WORK/chelpg_sampling_panel_identity_v1/collection_job_1201211.json" \
  --output "$ALQ_WORK/chelpg_sampling_review_report_v1"
```

## Preparation/execution interfaces used

The original qualification uses `scripts/mace_chelpg_sampling.py`: `prepare`,
`dry-run`, `execute`, `collect-native`, `prepare-identity`, `dry-run-identity`,
`execute-identity`. Its recorded configuration is
`chelpg_sampling_qualification_config_v1/config.json` under the workspace root.

The five-geometry extension uses `scripts/mace_chelpg_sampling_panel.py`:
`prepare --config ABSOLUTE_CONFIG --output NEW_DIRECTORY`, followed by
`dry-run --manifest ABSOLUTE_MANIFEST`. The actual fully specified configuration
is `chelpg_sampling_panel_config_v1/config.json`. The initial two property
results are reused only after exact source and density-identity checks.

For an explicitly prepared inventory, submit the pinned manifest to
`run_chelpg_sampling_panel.sbatch`; it executes four concurrent 16-rank native
replays. Then `prepare-identity --population ABSOLUTE_POPULATION --output
NEW_DIRECTORY` constructs actual density queries. Submit that manifest to
`run_chelpg_sampling_panel_identity.sbatch`. Both submit scripts require one
absolute manifest argument. No hidden input defaults or source edits are needed.

Do not resubmit the completed inventory merely to review it. `collect` preserves
partial/unavailable rows and starts no executable. Execution validates retained
NoIter receipts before reusing completed tasks; a partial/invalid attempt is
retained and requires an explicit fresh retry preparation. No failed score is
replaced with baseline data. This interface produces a numerical diagnostic,
not a calibrated biological classifier.
