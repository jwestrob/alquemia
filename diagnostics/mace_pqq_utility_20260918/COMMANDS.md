# PQQ utility benchmark operations

The timing_v2 jobs are already submitted: DFT1201562, MACE1201566. Do not launch
duplicates. Actual argv and scheduler responses are in the respective
`timing_v2/{DFT,MACE}_submission.json`. The failed MPI startup is preserved in
timing_v1; it is not successful timing evidence.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
PQQ_DRIVER=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
PQQ_WORK="$PWD/workspaces/mace_pqq_utility_20260918"
squeue -j 1201562,1201566
python diagnostics/mace_pqq_utility_20260918/test_benchmark.py
"$PQQ_DRIVER" diagnostics/mace_pqq_utility_20260918/test_source_selector.py
python diagnostics/mace_pqq_utility_20260918/benchmark.py collect \
  --manifest "$PQQ_WORK/timing_v2/manifest.json" \
  --output "$PQQ_WORK/timing_v2/comparison_v1.json"
```

Collect can report partial results but never fills missing scores or reports
partial-pair timings as the full campaign. Its output path must be new; retain
partial comparisons under new versioned names. Final reporting must include
Slurm allocation receipts and the `/usr/bin/time` files, so job-level startup
overhead and failed attempts are not lost.

After both jobs are terminal and all25pairs are complete, build the final
report (this refuses running/incomplete campaigns):

```bash
python diagnostics/mace_pqq_utility_20260918/report.py \
  --workspace "$PQQ_WORK" --output "$PQQ_WORK/final_report_v1"
```

The final report records literal frozen-band calls separately from reproduction
within the predeclared0.01tolerance, all100successful endpoint receipts, the
failed initial MPI startup, phase timing and full allocation costs.

The source-backed accuracy replay is already complete in accuracy_v1. To replay
without inference to a fresh location:

```bash
python diagnostics/mace_pqq_utility_20260918/accuracy.py \
  --root "$PWD" --output "$PQQ_WORK/accuracy_replay_v2"
```

Re-preparing a timing manifest alone launches nothing:

```bash
"$PQQ_DRIVER" diagnostics/mace_pqq_utility_20260918/benchmark.py prepare \
  --root "$PWD" --output "$PQQ_WORK/timing_reprepare_v3" \
  --frozen-mace "$PQQ_WORK/interface_source_v2/implementation"
```

Each scientific computation is delegated to the existing task runner. DFT
uses copied byte-identical templates/XYZ and allocation-specific `%pal`.
MACE uses the frozen two-call prepared-input interface with only the recorded
source-selector guard fix. No altered weights, adapters, coordinates or bands.
