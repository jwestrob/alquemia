# Saved responsive density / current GK functional: operations

Jobs 1201135 (eight saved-density utilities) and 1201137 (full model) completed.
All four ordering, the partition, and all 32 numerical checks pass. These are
two consumed development groups. See [result and limits](TRIAL_DENSITY_GK_REPORT.md)
and the pre-output [plan](TRIAL_DENSITY_GK_PLAN.md). The preceding failed candidate
is immutable. All paths below are explicit; dry-runs launch no scientific calls.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/trial_density_gk_fields_v1/implementation/mace_trial_density_fields.py dry-run --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/trial_density_gk_fields_v1/manifest.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/trial_density_gk_pipeline_v1/implementation/mace_density_gk_hybrid.py pipeline-dry-run --config /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/trial_density_gk_pipeline_v1/config.json
```

The dependency ran prepare, execute and collect through the existing hybrid
CLI after the real field result passes. It preserves the finite111-task logical
inventory, reuses37identical environment-only components, and executes34new
static energies,32field queries and40responses. It adds no DFT, MACE, fit,
geometry or cutoff change. Its frozen source/configuration lives in
`trial_density_gk_pipeline_v1`. The actual native executable/library is the same
qualified `density_gk_hybrid_software_v1` used by1201074.

Parent coordinates/key files are copied exactly, including actual rigid
geometries; the responsive source charge files are separate. Reuse compares
actual geometry, key, mask, charge, backend and native parameters. Response
reuse additionally demands identical supplied-field bytes. Historical native
costs remain attached to source receipts and are excluded from new-call sums.
No baseline result can fill a missing responsive correction.

## Actual completed preparation

- Source audit:8compatible saved quantum states; exact old-field subtraction
  reproduces archived intrinsic energies within1.165e-10kcal. Intrinsic
  polarization costs37.355376--41.323737kcal are allnonnegative. No new science.
- Responsive boundary:16actual native initializations and12paired reference/
  cavity checks pass. Same exterior environments as the previous candidate.
  Preparation65.835947813s wall/47.364624084CPU; native16.733242134wall/16.553698CPU.
- Query preparation:13.711318422s wall/12.687261401CPU,49observations per actual
  environment site. No new potential value was manufactured or inferred.
- Existing field4tests pass4.758s, boundary3pass39.440s, hybrid6pass32.976s.
  All 12 distinct new tests now pass, including the actual utility and full-model
  integration tests. No final test remains skipped. Five source/boundary/reuse/
  pipeline tests pass across focused runs. One erroneous inventory assertion(29instead of37)
  was corrected from the actual task list before any new energy/response call.

## Check and recollect completed actual outputs

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/tests -p test_mace_trial_density.py -v
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/tests -p test_mace_trial_density_fields.py -v
```

Completed actual results:
`trial_density_gk_fields_v1/report_job_1201135/result.json` and
`trial_density_gk_hybrid_v1/collection_job_1201137.json`.
The pipeline phase/failure receipt is
`trial_density_gk_pipeline_v1/pipeline_receipt_job_1201137.json`.
Read the Slurm logs if the prerequisite/prepare phase fails; no numerical score
is supplied in that state. Collection keeps all original identity/refinement/
radius/partition/ordering tolerances. A trial density is not a self-consistent
AMOEBA/GK quantum solution, and the new score has no absolute calibration.

For read-only recollection after a real manifest exists, write to a fresh path:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/trial_density_gk_hybrid_v1/implementation/mace_density_gk_hybrid.py collect --manifest /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/trial_density_gk_hybrid_v1/manifest.json --output /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_omol_20260917/trial_density_gk_hybrid_v1/collection_replay_01.json
```

For recovery only after confirming there is no live executor, the queued
pipeline can be resubmitted with the same frozen config and recorded resource
wrapper. `execute --retry-failed` is available on an existing prepared manifest
for explicit technical recovery with unchanged science. Do not alter source
states or overwrite preparation/results to make a retry succeed.

## Compare the matched archived candidates

This reads actual results and writes a fresh comparison; it launches no solver.
From the repository root:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_density_gk_compare.py --parent workspaces/mace_omol_20260917/density_gk_hybrid_v1/collection_job_1201074.json --trial workspaces/mace_omol_20260917/trial_density_gk_hybrid_v1/collection_job_1201137.json --output workspaces/mace_omol_20260917/trial_density_gk_comparison_review.json
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_mace_trial_density_gk.py -v
```

Do not re-execute completed jobs to regenerate a report. For a genuinely new
experiment, create a new versioned configuration/manifest with explicit sources
and scope; do not edit the frozen trial config or overwrite its scientific output.
