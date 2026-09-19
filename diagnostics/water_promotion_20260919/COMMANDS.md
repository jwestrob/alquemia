# Promoted baseline operations

Run from the repository root. The default water policy is `contextual_if_supported`.
Use `--water-policy original` at request creation to select the historical inputs.
The release is pinned in `params/baseline_water_v1.json`; no source edits are needed.

## Reproduce the completed release checks, without new calculations

```bash
BW_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$BW_PY" scripts/affordable_workflow.py baseline dry-run \
  --manifest workspaces/water_promotion_20260919/dft_v2/manifest.json
"$BW_PY" scripts/affordable_workflow.py baseline dry-run \
  --manifest workspaces/water_promotion_20260919/dry_dft_v3/manifest.json
"$BW_PY" -m unittest discover -s tests -p test_baseline_water.py -v
```

## Request, prepare, execute, collect and report

This concrete example reuses the existing alpha-lactalbumin demonstration. All
output paths are new and write-once; it does not rerun its successful chemistry.
For a new sample supply its supported whole-protein `preparation.json` through
`--preparation`, omit `--reuse-request`, and give new output paths. An input
requires the existing source-backed amide-v3 core and explicit assembly mapping;
a raw CIF alone is not this interface. Dry canonical PQQ is an exact identity.
Use `--context-group` only for explicitly matched, same-indexed constructs;
independent samples otherwise receive independent contexts.

```bash
"$BW_PY" scripts/affordable_workflow.py baseline request \
  --preparation workspaces/mace_global_benchmark_20260916/prepared_v1/ALPHA_1F6S/preparation.json \
  --preparation workspaces/mace_global_benchmark_20260916/prepared_v1/ALPHA_6IP9/preparation.json \
  --context-group ALPHA_1F6S=alpha_same_indexed_construct \
  --context-group ALPHA_6IP9=alpha_same_indexed_construct \
  --reuse-request diagnostics/contextual_water_20260919/REPLAY_REQUEST.json \
  --output workspaces/water_promotion_20260919/operator_request_v1.json
"$BW_PY" scripts/affordable_workflow.py baseline prepare \
  --request workspaces/water_promotion_20260919/operator_request_v1.json \
  --output workspaces/water_promotion_20260919/operator_prepared_v1
"$BW_PY" scripts/affordable_workflow.py baseline water-execute \
  --plan workspaces/water_promotion_20260919/operator_prepared_v1/plan.json
"$BW_PY" scripts/affordable_workflow.py baseline prepare-dft \
  --plan workspaces/water_promotion_20260919/operator_prepared_v1/plan.json \
  --reuse-manifest workspaces/benchmark_set_20260915/ready_tasks_v4/manifest.json \
  --reuse-manifest workspaces/hydration_network_20260918/core_transfer_v1/manifest.json \
  --output workspaces/water_promotion_20260919/operator_dft_v1
"$BW_PY" scripts/affordable_workflow.py baseline dry-run \
  --manifest workspaces/water_promotion_20260919/operator_dft_v1/manifest.json
"$BW_PY" scripts/affordable_workflow.py baseline execute \
  --manifest workspaces/water_promotion_20260919/operator_dft_v1/manifest.json
"$BW_PY" scripts/affordable_workflow.py baseline collect \
  --manifest workspaces/water_promotion_20260919/operator_dft_v1/manifest.json \
  --output workspaces/water_promotion_20260919/operator_result_v1.json
"$BW_PY" scripts/affordable_workflow.py baseline report \
  --collection workspaces/water_promotion_20260919/operator_result_v1.json \
  --output workspaces/water_promotion_20260919/operator_report_v1.md
```

For actual new tasks, `water-execute` belongs in a suitable GPU allocation using
the existing pinned MACE runner; `execute` belongs in a CPU Slurm allocation
with at least64 CPUs (four concurrent16-rank ORCA workers, fewer workers when
fewer tasks remain). Add the installed ORCA directory to PATH as required by
the existing runner. These commands do not submit jobs or choose a partition.
The finite rendered manifest exposes all work before execution. There is no
project runtime budget stopping rule. In the replay above both execute operations
are zero-work and require no allocation.

Each wet case retains both original and prepared endpoint pairs. Existing original
outputs can be reused via repeated `--reuse-manifest`; without compatible archived
outputs both pairs need calculation. This cost is explicit, not hidden in preparation.
Dry cases need only the unchanged original pair. Unsupported water-bearing cores,
missing nearby chemistry and unconverged proposals fail explicitly.

A partially executed DFT manifest can be resumed with `execute` if completed tasks
are intact and missing tasks have no partial files. Preserve failed/partial outputs
and create a fresh retry with verified successful results reused:

```bash
"$BW_PY" scripts/affordable_workflow.py baseline prepare-retry \
  --manifest workspaces/water_promotion_20260919/operator_dft_v1/manifest.json \
  --output workspaces/water_promotion_20260919/operator_retry_v1
```

Then use `dry-run`, `execute`, `collect` and `report` on the returned manifest.
A missing prepared score remains unavailable even when its original score exists.
Changed wet scores have no inherited PQQ bands or aquo reference. No entropy,
occupancy probability, water insertion or per-case favorable-score selection occurs.
