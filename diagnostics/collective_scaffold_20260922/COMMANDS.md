# Collective scaffold mechanics operations

The scientific scope is the fixed common eight sources and three donor targets
per source. Preserve the failed CUDA records in `searches_v1` and the superseded
achiral-guard attempts in `searches_v2`; only corrected `searches_v3` is eligible
for scoring. Do not resubmit a completed search or use this
mechanics branch to run MACE/GFN2 scoring.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
SCAFFOLD_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
SCAFFOLD_RUN=$PWD/workspaces/collective_scaffold_20260922/searches_v3
"$SCAFFOLD_PY" "$SCAFFOLD_RUN/implementation/scaffold_accommodation.py" validate --manifest "$SCAFFOLD_RUN/manifest.json"
"$SCAFFOLD_PY" -m unittest discover -s tests -p test_scaffold_accommodation.py -v
```

Preparation accepts explicit paths and writes a new immutable directory. This is
geometry-only; the destination below must not already exist:

```bash
"$SCAFFOLD_PY" scripts/scaffold_accommodation.py prepare \
  --parents workspaces/collective_scaffold_20260922/parents_v1/manifest.json \
  --inputs diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json \
  --agreement diagnostics/scaffold_restart_round_20260922/PLAN.md \
  --target-mapping diagnostics/scaffold_restart_round_20260922/TARGET_MAPPING.md \
  --plan diagnostics/collective_scaffold_20260922/EXECUTION_PLAN.md \
  --technical-addendum diagnostics/collective_scaffold_20260922/STEREOCHEMISTRY_RECOVERY.md \
  --output workspaces/collective_scaffold_20260922/preparation_replay
```

Actual submission is recorded in the run's `SUBMISSION.json`, with the exact
wrapper and manifest hashes. `run_mechanics.sbatch` executes only prepared tasks,
preserves existing receipts and incomplete attempts, and writes a terminal
collection automatically. Its allocation is one H200, 32 CPUs and 200000 MiB.

To collect current or final existing outputs into a fresh artifact:

```bash
"$SCAFFOLD_PY" "$SCAFFOLD_RUN/implementation/scaffold_accommodation.py" collect \
  --manifest "$SCAFFOLD_RUN/manifest.json" \
  --output "$SCAFFOLD_RUN/collection_replay.json"
```

The collection preserves all 24 declared targets. Each complete receipt pins full
parent coordinates, regenerated Ca/La contexts, actual initial/final energies,
every force evaluation, physical checks and stationarity status. Missing or
failed targets are unavailable. A feasible finite proposal is not automatically
a stationary protein structure. Parent FF work is never a score contribution.

After the real terminal collection exists, reproduce the geometry/constraint
summary without any force call, using a fresh output destination:

```bash
"$SCAFFOLD_PY" diagnostics/collective_scaffold_20260922/summarize_mechanics.py \
  --collection "$SCAFFOLD_RUN/collection_1210206.json" \
  --output "$SCAFFOLD_RUN/mechanics_summary_replay.json"
```
