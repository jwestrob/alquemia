# Fast PQQ standard mode

Run from the repository root on biotite. These examples use real source inputs;
no prior MACE energies or user edits to source code are required. The explicit
source request specifies assembly, metal/PQQ selectors and the five canonical
homologous residue roles. This is not an automatic raw metagenome/fold annotation
service. A structurally compatible unseen source retains `unvalidated_input_domain`.

```bash
export PQQ_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Prepare a source request

Use a new output directory for each execution. The existing `standard_1H4I_v1`
directory contains the release check; the command below uses a fresh example
output and is not intended to be rerun against an existing directory.

```bash
"$PQQ_PY" scripts/affordable_workflow.py standard prepare \
  --request diagnostics/pqq_fast_release_20260920/examples/1H4I_source.json \
  --output workspaces/pqq_fast_release_20260920/example_1H4I
"$PQQ_PY" scripts/affordable_workflow.py standard validate \
  --plan workspaces/pqq_fast_release_20260920/example_1H4I/plan.json
```

Default `--mode auto` selects `fast_PQQ_OMOL_GFN2_ALPB_v1` for a supported explicit
PQQ source request. `--mode fast-pqq` requires that route. Both reject unsupported
source chemistry and retain failed component status; neither substitutes a DFT
or old cached score after failure.

## Execute the finite source-to-score plan

```bash
sbatch \
  --output=workspaces/pqq_fast_release_20260920/example_1H4I/slurm_%j.out \
  --error=workspaces/pqq_fast_release_20260920/example_1H4I/slurm_%j.err \
  diagnostics/pqq_fast_release_20260920/run_standard.sbatch \
  "$PWD/workspaces/pqq_fast_release_20260920/example_1H4I/plan.json"
```

This uses one H200, 32 CPUs/MPI slots and 200000 MiB. It runs fresh source
preparation, two native MACE energies and four native GFN2 energies per site,
then collects `result_JOBID.json`. Existing runner snapshots, model and source
pins are verified. The job manifest and actual receipts are retained. No
folding, full-protein simulation, geometry search or DFT runs are included.

## Keep the DFT reference explicitly available

Use the same supported source request with `--mode dft-reference` and a fresh
output directory. It emits the exact original canonical-core r2SCAN-3c/CPCM
recipe and original DFT reference/bands. The standard entrypoint uses the same
batch wrapper for two 16-rank native DFT endpoints; this release did not rerun
those DFT energies. This branch's source request must still meet the supported
complete-context preparation policy. Legacy prepared-core DFT remains available
through the unchanged `baseline` entrypoint.

```bash
"$PQQ_PY" scripts/affordable_workflow.py standard prepare \
  --request diagnostics/pqq_fast_release_20260920/examples/1H4I_source.json \
  --mode dft-reference \
  --output workspaces/pqq_fast_release_20260920/example_1H4I_DFT
```

## Other chemistry

Existing `baseline_contextual_water_v1` requests sent to `standard prepare`
produce the unchanged staged water-preparation/DFT plan and identify the next
stage. Continue using `affordable_workflow.py baseline` for its existing resource
stages and collection. `standard execute` reports that staged entrypoint; it
does not pretend to run those stages or silently reinterpret a generic site as
PQQ. Existing baseline operations, raw inbox watchers and PLM results are intact.

## Inspect the completed release without computing

```bash
"$PQQ_PY" scripts/affordable_workflow.py standard validate \
  --plan workspaces/pqq_fast_release_20260920/standard_1H4I_v1/plan.json
"$PQQ_PY" -c 'import json; r=json.load(open("diagnostics/pqq_fast_release_20260920/RESULT.json")); print({k:r[k] for k in ("promoted","calibration_correct","transfer_correct","allocated_core_seconds","allocated_GPU_seconds")})'
```

The standard CLI result is
`workspaces/pqq_fast_release_20260920/standard_1H4I_v1/result_1203729.json`.
The full 28-case source/preparation/scoring collection is pinned by
`diagnostics/pqq_fast_release_20260920/RELEASE_PANEL.json`. These are completed
runs. Read/validate them rather than submit their manifests again.

To collect a newly completed standard run, provide its plan and a new output
path to `standard collect`; the batch wrapper already does that automatically.
Existing files are immutable. For a failed attempt, retain its directory and
receipts and use a fresh explicit output directory for a same-input recovery.
