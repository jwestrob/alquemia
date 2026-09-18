# Resolved-boundary component operations

These are opt-in development operations. They do not run a production scorer.
Existing outputs are immutable; execution refuses to overwrite an attempt.
All paths below refer to the declared real GGR inventory.

```bash
DDX_BASE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
DDX_PROJECT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
DDX_WORK_ROOT="$DDX_PROJECT/workspaces/mace_omol_20260917/ddx_source_recovery_v1"
export OPENBLAS_NUM_THREADS=1
```

Audit / dry-run, without importing or calling the native solver:

```bash
"$DDX_BASE_PY" "$DDX_WORK_ROOT/implementation/mace_ddx_recovery.py" dry-run \
  --manifest "$DDX_WORK_ROOT/manifest.json"
```

Preparation was executed with the recorded configuration:

```bash
"$DDX_BASE_PY" "$DDX_PROJECT/scripts/mace_ddx_recovery.py" prepare \
  --config "$DDX_PROJECT/workspaces/mace_omol_20260917/ddx_source_recovery_config_v1/config.json" \
  --output "$DDX_WORK_ROOT"
```

It pins18missing solve roles and2reuses across8groups. The fully reusable coarse
2FW0 group was collected locally with zero native calls. Seven group jobs
1201286–1201292 were submitted using `run_ddx_recovery.sbatch`; exact commands,
wrapper hash and job receipts are in `submissions.json` and the per-group
submission files. Do not resubmit those live jobs. The wrapper accepts the
absolute manifest and an exact group ID; each solver uses64CPUs/128GiB with
maxiter1200 and no automatic retry. There is no project time/compute budget.

Next operation, after all seven group receipts exist:

```bash
"$DDX_BASE_PY" "$DDX_WORK_ROOT/implementation/mace_ddx_recovery.py" collect \
  --manifest "$DDX_WORK_ROOT/manifest.json" \
  --output "$DDX_WORK_ROOT/collection.json"
"$DDX_BASE_PY" "$DDX_WORK_ROOT/implementation/mace_ddx_source_report.py" \
  --collection "$DDX_WORK_ROOT/collection.json" \
  --output "$DDX_PROJECT/workspaces/mace_omol_20260917/ddx_source_recovery_report_v1"
```

Collection/reporting starts no solver. It checks actual coefficients, units,
source identity, reuse identity, runtime parameters and the frozen physical
screens. Missing/failed energies remain unavailable. Ca-minus-La and matched
GGR structural differences are solvent-component contrasts, not calibrated
discriminator scores. The coarsest2FW0 reciprocity check already fails0.05kcal
(observed0.122802); retain it even if finer representations improve.

Pinned real-fixture tests (no fabricated solver output):

```bash
cd "$DDX_PROJECT"
"$DDX_BASE_PY" -m unittest discover -s tests -p 'test_mace_ddx_*.py' -v
```

The original incomplete20-role run and the successful single-state logged
recovery remain under `ddx_source_self_v2` and `ddx_convergence_v1`. Their
records are not changed by the new collection. See DDX_SOURCE_SELF_PLAN,
DDX_CONVERGENCE_REPORT and DDX_ITERATION_RECOVERY_PLAN for the actual scope.


## Current resolution and polarization diagnostics

These commands collect completed existing work; they do not submit new jobs.
Run from the repository root. All scientific configurations are in the pinned
manifests and their implementation snapshots.

```bash
OPENBLAS_NUM_THREADS=1 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/mace_omol_20260917/ddx_resolution_v1/implementation/mace_ddx_resolution.py \
  dry-run --manifest workspaces/mace_omol_20260917/ddx_resolution_v1/manifest.json

# Job 1201312 writes collection_job_1201312.json on completion. For a separate
# read-only collection after all four groups finish, use a new output filename:
OPENBLAS_NUM_THREADS=1 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/mace_omol_20260917/ddx_resolution_v1/implementation/mace_ddx_resolution.py \
  collect --manifest workspaces/mace_omol_20260917/ddx_resolution_v1/manifest.json \
  --output workspaces/mace_omol_20260917/ddx_resolution_v1/recollection.json

OPENBLAS_NUM_THREADS=1 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/mace_omol_20260917/frozen_response_v2/implementation/mace_frozen_response.py \
  dry-run --manifest workspaces/mace_omol_20260917/frozen_response_v2/manifest.json
```

Corrected native-functional results are in
`frozen_response_v2/collection_job_1201310.json`. They pass their accounting
checks; no conductor transfer or biological score is present. The earlier
wrapper failure in `frozen_response_v1` is retained. Completed PCM/conductor
inventories and their numerical reports are described in CURRENT.md.
