# Experimental three-source command, interface v2

The owned entrypoint remains `scripts/pqq_three_source_envelope_execution.py`.
New `prepare` calls create interface-v2 plans from one or more explicitly named,
complete prepared triples. Source IDs and protein names are not a whitelist.
Model/geometry/solvent policy and the frozen envelope reference are unchanged.
No command below was submitted for molecular execution in this API task.

## Prepare and inspect a real single-protein request

The source request/preparation route is documented in
[the preparation commands](../pqq_three_source_envelope_20260923/COMMANDS.md).
It accepts explicit source/config/reference paths and currently requires exact
compatible archived protonation. It checks all three structures' actual protein
and physical state; the physical anchor does not require a historical label.

This runnable example creates a new finite plan for the already prepared07ab
triple. All paths are explicit. Use a new output directory for each immutable plan.

```bash
py=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$py" scripts/pqq_three_source_envelope_execution.py prepare \
  --preparations workspaces/pqq_three_source_envelope_20260923/plm_v1/PQQSEQ_07ab500e3df76b30d71c/PREPARATION.json \
  --strict-qualification workspaces/strict_native_pool_20260923/run_v1/COMPARISON.json \
  --rank-qualification workspaces/native_gfn2_rank_panel_20260923/run_v2/COMPARISON.json \
  --maxiter-qualification workspaces/adaptive_maxiter_20260922/prepared_v1/collection_1209876.json \
  --agreement diagnostics/pqq_three_source_envelope_api_20260923/PLAN.md \
  --output workspaces/pqq_three_source_envelope_api_20260923/command_example_v1
"$py" scripts/pqq_three_source_envelope_execution.py dry-run \
  --plan workspaces/pqq_three_source_envelope_api_20260923/command_example_v1/plan.json
```

For both explicit groups, pass both actual `PREPARATION.json` paths after
`--preparations`; the input list determines the finite plan. Each group must be
complete. Unsupported or mismatched groups cause preparation to fail before
staging new molecular inputs. No member is silently dropped.

Immediate read-only inspection of a completed v2 preflight:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/pqq_three_source_envelope_execution.py dry-run \
  --plan workspaces/pqq_three_source_envelope_api_20260923/preflight_v1/07ab/plan.json
```

## Execution operation (available, not run in this task)

The existing pinned wrapper can execute a declared plan in its required
32CPU/oneH200/200000MiB allocation. This is an explicit single-protein example,
not a default or cohort submission. Actual source/reference eligibility still
applies; scientific use remains experimental.

```bash
root=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
plan="$root/workspaces/pqq_three_source_envelope_api_20260923/preflight_v1/07ab/plan.json"
sbatch --parsable --partition=gpu --nodelist=node-224-2t-8gpu-1 \
  --nodes=1 --ntasks=32 --cpus-per-task=1 --gres=gpu:1 --mem=200000M \
  --job-name=nikasha-envelope07 \
  --output="$root/workspaces/pqq_three_source_envelope_api_20260923/preflight_v1/07ab/slurm_%j.out" \
  --error="$root/workspaces/pqq_three_source_envelope_api_20260923/preflight_v1/07ab/slurm_%j.err" \
  "$root/diagnostics/pqq_three_source_envelope_execution_20260923/run.sbatch" "$plan"
```

The wrapper uses `execute --plan`, then always attempts `collect` and `report`,
retaining the molecular exit status. It does not silently retry an already
started plan. Results retain all source members, original/accommodated contrasts,
both common-pool selections/works, strict median/spread and missing reasons.

## Historical result replay

The completed six-source result remains v1. Current validation and read-only
collection support it, while its immutable executed script remains unchanged.
There is no new v2 molecular result or implied v2 rerun.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/pqq_three_source_envelope_execution.py report \
  --result workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RESULT_1211626.json \
  --output workspaces/pqq_three_source_envelope_api_20260923/HISTORICAL_v1_REPORT.md
```

Per complete triple:6 fresh MACE origins,6 bounded searches,≤6crossMACE and36
strict scalar cells. Search-internal evaluation count is data dependent. These
counts multiply by the explicitly declared group count; they do not allocate or
launch work. The operative SLSQP tolerance remains `optimizer_ftol=1e-8` on the
Hartree-equivalent scaled objective, with multiplier0.03674932217934773; inherited
`ftol_eV=1e-9` metadata is retained byte-for-byte. Do not describe the operative
criterion as an unscaled1e-8eV tolerance.
