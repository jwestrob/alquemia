# Reciprocal metal-position MACE scoring

Prepared implementation; no model loading, molecular call or scientific test has
run during this code preparation. Bind only the actual Slurm preparation result.

`SCORE.py` calls the unchanged pinned `CPU_BATCH.WarmCPU` and `cell` implementation.
One persistent CPU calculator performs exactly eight transplanted evaluations.
Its ordinary energy call also caches forces; there is one native model forward
per endpoint, no optimizer and no additional force evaluation. All eight original
192-atom results are reused through the pinned ExaA donor-response inventory;
the collector checks request XYZ, charge, multiplicity, model, precision, runtime
versions and actual execution evidence before comparing them.

The checkpoint is `MACE-omol-0-extra-large-1024.model`, SHA256
`9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`.
The original OMOL head, float64, charge −3/−2 and singlet multiplicity inputs are
unchanged. Packing job4091 already qualified this seven-thread native CPU path
against both archived Ca/La force/energy fixtures. Its actual qualification and
worker summary are pinned and reused; no fixture calls are added.

## Root commands

Metadata-only binding after Curie's actual preparation completes:

```bash
/usr/bin/python3 diagnostics/exaa_metal_position_20261014/scoring/SCORE.py prepare \
  --preparation /ABS/PREPARATION.json --output /ABS/FRESH_SCORE_WORKSPACE
```

This prints the actual `CONFIG.json` pin and snapshots this adapter and launcher.
Root substitutes those exact paths and hash below:

```bash
sbatch --parsable --nodelist=node-64-768g-5 \
  /ABS/FRESH_SCORE_WORKSPACE/RUN.sbatch \
  /ABS/FRESH_SCORE_WORKSPACE/CONFIG.json CONFIG_SHA run
```

Use `preflight` instead of `run` for a separate seven-CPU Slurm preflight without
loading the model. The normal `run` performs that same preflight before the eight
calls, so a separate preflight is optional. Execution is immutable: an existing
execution directory refuses a second run. Completed cell receipts and partial
collection remain available after failures; no hidden retry or fallback occurs.

## Allocation and measured precedent

The launcher requests one node/seven tasks, without memory/time overrides. The
current shared default gives about28GB; packing4091 measured8,931,864KiB peak
hostRSS. Its fourteen native forwards on larger249-atom geometries consumed
243.43s worker wall time including1.90s model loading. These eight192-atom calls
are expected to take roughly2–3minutes plus preflight; this is an estimate from
actual receipts, not measured throughput for the new task.

The parent Python is `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python`.
It keeps numerical thread variables at1 while the persistent model subprocess
uses the existing pinned `workspaces/mace_hybrid_20260916/software_v1/venv/bin/python`
with seven Torch threads and one inter-op thread. The parent waits during model
calls. No model/environment installation or production edit is required.

`execution/FINAL_COLLECTION.json` retains both original/transplanted energies,
per-metal origin→adaptive_La work, Ca-minus-La differential work, each host's
delta_D, and the target-minus-control contrast. Missing cells remain explicit.
No solvent term, classifier band, optimization or biological label is inferred.
