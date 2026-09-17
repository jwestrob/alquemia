# Charge/dipole source operations

V1 completed as1201022. The denser-sampling V2 has a separate manifest and
submission receipt under `distributed_source_v2_recovery_v1`. Check that receipt
and the live queue before executing anything; do not duplicate a running job.
The initial `distributed_source_v2` preparation failed before any fit/native
call due to a local variable shadowing a function. It is preserved with its
failure record and original implementation; the recovery changes no physics.

```bash
ALQUEMIA_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ALQUEMIA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
SOURCE_BASE=$ALQUEMIA_ROOT/workspaces/mace_omol_20260917
SOURCE_V1=$SOURCE_BASE/distributed_source_v1
SOURCE_V2=$SOURCE_BASE/distributed_source_v2_recovery_v1

"$ALQUEMIA_PY" "$SOURCE_V1/implementation/mace_distributed_source.py" dry-run --manifest "$SOURCE_V1/manifest.json"
"$ALQUEMIA_PY" "$SOURCE_V2/implementation/mace_distributed_source.py" dry-run --manifest "$SOURCE_V2/manifest.json"
# Read-only science recollection; new report output is exclusive.
"$ALQUEMIA_PY" "$SOURCE_V2/implementation/mace_distributed_source.py" collect \
 --manifest "$SOURCE_V2/manifest.json" --output "$SOURCE_V2/recollection_v1"
cd "$ALQUEMIA_ROOT"
"$ALQUEMIA_PY" -m unittest discover -s tests -p test_mace_distributed_source.py -v
```

An input replay, if independently needed, performs no fits or native calls:

```bash
"$ALQUEMIA_PY" "$ALQUEMIA_ROOT/scripts/mace_distributed_source.py" prepare \
 --fields "$SOURCE_BASE/qm_electric_field_v1/report_job_1201017/result.json" \
 --center-manifest "$SOURCE_BASE/explicit_field_short_v1/manifest.json" \
 --additional-training-report "$SOURCE_V1/report_job_1201022/result.json" \
 --plan "$ALQUEMIA_ROOT/diagnostics/mace_omol_20260917/DISTRIBUTED_SOURCE_SAMPLING_PLAN.md" \
 --output "$SOURCE_BASE/distributed_source_v2_replay_v1"
# Optional explicit execution of that new replay:
# sbatch "$ALQUEMIA_ROOT/diagnostics/mace_omol_20260917/run_distributed_source.sbatch" \
#  "$SOURCE_BASE/distributed_source_v2_replay_v1/manifest.json"
```

Omitting `--additional-training-report` selects V1; use DISTRIBUTED_SOURCE_PLAN.md
as its plan and a new V1 output directory. Each version retains its own source,
sampling, fit and implementation hashes. All successful identical tasks are
reused on restart. Failed attempts remain visible and need explicit
`execute --retry-failed` in the declared allocation; any repeated fit/native
call counts in the execution history. No reported scalar is an affinity score.
