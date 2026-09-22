# Three-point native GFN2 numerical diagnostic

From the repository root, with the existing CPU environment:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
DIAG_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
DIAG_ROOT="$PWD/workspaces/adaptive_maxiter_20260922/prepared_v1"
"$DIAG_PY" "$DIAG_ROOT/implementation/adaptive_maxiter_diagnostic.py" validate --manifest "$DIAG_ROOT/manifest.json"
```

Preparation already executed, retaining all three exact source tasks:

```bash
"$DIAG_PY" scripts/adaptive_maxiter_diagnostic.py prepare \
  --failed-manifest workspaces/adaptive_accommodation_20260922/common_pool_v1/solvent/shard_0/manifest.json \
  --pilot-manifest workspaces/nikasha_shared_pool_20260922/pilot_v2/manifest.json \
  --agreement diagnostics/adaptive_maxiter_20260922/PLAN.md \
  --output workspaces/adaptive_maxiter_20260922/prepared_v1
```

The output directory is write-once. Job **1209876** already executes the exact
manifest below. Do not resubmit it. Its wrapper automatically collects terminal
successes/failures into `collection_1209876.json`.

```bash
sbatch --parsable \
  --output="$DIAG_ROOT/slurm_%j.out" \
  --error="$DIAG_ROOT/slurm_%j.err" \
  diagnostics/adaptive_maxiter_20260922/run.sbatch "$DIAG_ROOT/manifest.json"
```

To collect a new read-only snapshot, choose a new output filename:

```bash
"$DIAG_PY" "$DIAG_ROOT/implementation/adaptive_maxiter_diagnostic.py" collect \
  --manifest "$DIAG_ROOT/manifest.json" \
  --output "$DIAG_ROOT/collection_review_v1.json"
```

Manifest SHA-256:
`018f68734e0817020e0e78076752c5d4434da359a5401e9e0b85c4b1028a4809`.
This diagnostic creates no recovery overlay or score. Its original failed pool
cell remains unavailable in the primary record.
