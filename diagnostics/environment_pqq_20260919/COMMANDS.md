# Compact-context PQQ operations

Run from the repository root. The 52-endpoint science job is complete. These
commands expose the original operations; do not rerun completed science as a
status check. Existing outputs are immutable and creation refuses collisions.

```bash
ENV_PQQ_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
ENV_PQQ_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ENV_PQQ_ROOT"
```

Validate the actual frozen preparation and reuse without calculations:

```bash
"$ENV_PQQ_PY" scripts/environment_pqq_context.py validate --manifest "$ENV_PQQ_ROOT/workspaces/environment_pqq_20260919/prepared_v1/manifest.json"
"$ENV_PQQ_PY" -m unittest discover -s tests -p test_environment_context_chemistry.py -v
```

Reproduce the inventory and preparation in new paths, with no energy calls:

```bash
"$ENV_PQQ_PY" scripts/environment_pqq_context.py inventory --source "$ENV_PQQ_ROOT/diagnostics/mace_omol_20260917/result.json" --static-config "$ENV_PQQ_ROOT/diagnostics/second_shell_20260919/CONFIG.json" --chemistry-version v2 --output "$ENV_PQQ_ROOT/workspaces/environment_pqq_20260919/inventory_replay_v1.json"
"$ENV_PQQ_PY" scripts/environment_pqq_context.py prepare --inventory "$ENV_PQQ_ROOT/workspaces/environment_pqq_20260919/inventory_replay_v1.json" --agreement "$ENV_PQQ_ROOT/diagnostics/environment_pqq_20260919/PLAN.md" --output "$ENV_PQQ_ROOT/workspaces/environment_pqq_20260919/preparation_replay_v1"
```

Collect and report existing actual outputs into fresh files (no new science):

```bash
"$ENV_PQQ_PY" scripts/environment_pqq_context.py collect --collection "$ENV_PQQ_ROOT/workspaces/environment_pqq_20260919/prepared_v1/collection_1202474.json" --output "$ENV_PQQ_ROOT/workspaces/environment_pqq_20260919/result_replay_v1.json"
"$ENV_PQQ_PY" scripts/environment_pqq_report.py --result "$ENV_PQQ_ROOT/workspaces/environment_pqq_20260919/result_replay_v1.json" --output "$ENV_PQQ_ROOT/workspaces/environment_pqq_20260919/report_replay_v1.json"
```

The exact finite execution interface is the pinned `run_mace.sbatch` and
`SUBMISSION.json`. The original submission was equivalent to this command;
**it has already run as1202474** and is recorded for reproducibility:

```bash
sbatch --output "$ENV_PQQ_ROOT/diagnostics/environment_pqq_20260919/mace_%j.out" --error "$ENV_PQQ_ROOT/diagnostics/environment_pqq_20260919/mace_%j.err" "$ENV_PQQ_ROOT/diagnostics/environment_pqq_20260919/run_mace.sbatch" "$ENV_PQQ_ROOT/workspaces/environment_pqq_20260919/prepared_v1/manifest.json"
```

A completed execution directory cannot be silently reused. New scientific
inputs need a new manifest/preparation version and new model-specific bands.
No command changes the default scorer or performs production rescoring.
