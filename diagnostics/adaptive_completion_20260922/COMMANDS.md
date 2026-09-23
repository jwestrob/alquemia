# Scaled angular proposal operations

Commands run from this repository. All scientific source paths are explicit.
Original30 jobs1209963/1209968 are complete; their saved results should be reused.
The primary225 manifest is prepared and **unlaunched**. Its execution command
below is for the coordinated next phase only, after root freezes the reference
and decides whether transfer is useful. Do not submit it automatically.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ADAPTIVE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
ADAPTIVE_ROOT=$PWD/workspaces/adaptive_completion_20260922
```

Read-only audit/dry-run using the exact prepared implementation snapshots:

```bash
"$ADAPTIVE_PY" "$ADAPTIVE_ROOT/original30_v1/implementation/adaptive_completion.py" dry-run --manifest "$ADAPTIVE_ROOT/original30_v1/manifest.json"
"$ADAPTIVE_PY" "$ADAPTIVE_ROOT/primary225_v1/implementation/adaptive_completion_folds.py" dry-run --manifest "$ADAPTIVE_ROOT/primary225_v1/manifest.json"
```

Existing collections:

```bash
cat "$ADAPTIVE_ROOT/original30_v1/TERMINAL_SUMMARY.json"
cat "$ADAPTIVE_ROOT/original30_v1/ACCOUNTING.txt"
cat "$ADAPTIVE_ROOT/primary225_v1/PREFLIGHT.json"
```

Collect a new snapshot without molecular calls; each destination must be new:

```bash
ADAPTIVE_SNAPSHOT=$(date -u +%Y%m%dT%H%M%SZ)
"$ADAPTIVE_PY" "$ADAPTIVE_ROOT/original30_v1/implementation/adaptive_completion.py" collect --manifest "$ADAPTIVE_ROOT/original30_v1/manifest.json" --output "$ADAPTIVE_ROOT/original30_v1/collection_$ADAPTIVE_SNAPSHOT.json"
"$ADAPTIVE_PY" "$ADAPTIVE_ROOT/primary225_v1/implementation/adaptive_completion_folds.py" collect --manifest "$ADAPTIVE_ROOT/primary225_v1/manifest.json" --output "$ADAPTIVE_ROOT/primary225_v1/collection_$ADAPTIVE_SNAPSHOT.json"
```

Exact preparation inputs, if a separate recorded recreation is required (do not
replace the existing manifests or submit this recreation as a new experiment):

```bash
ADAPTIVE_REPLAY=$(date -u +%Y%m%dT%H%M%SZ)
"$ADAPTIVE_PY" scripts/adaptive_completion.py prepare --angular-manifest "$PWD/workspaces/adaptive_accommodation_20260922/proposals_v1/manifest.json" --angular-manifest "$PWD/workspaces/adaptive_accommodation_20260922/remaining26_proposals_v1/manifest.json" --agreement "$PWD/diagnostics/adaptive_completion_20260922/PLAN.md" --output "$ADAPTIVE_ROOT/replay_original30_$ADAPTIVE_REPLAY"
"$ADAPTIVE_PY" scripts/adaptive_completion_folds.py prepare --source "$PWD/workspaces/accommodation_nonlinear_20260920/fold_proposals_v1/manifest.json" --base-collection "$PWD/workspaces/nikasha_shared_pool_20260922/primary225_v1/final_collection.json" --diagnostic "$PWD/workspaces/adaptive_accommodation_20260922/force_projection_v1.json" --agreement "$PWD/diagnostics/adaptive_completion_20260922/TRANSFER_PREPARATION.md" --output "$ADAPTIVE_ROOT/replay_primary225_$ADAPTIVE_REPLAY"
```

Future primary225 execution inside the agreed one-H200,32-CPU,200000-MiB Slurm
allocation, once root explicitly releases this already prepared phase:

```bash
"$ADAPTIVE_PY" "$ADAPTIVE_ROOT/primary225_v1/implementation/adaptive_completion_folds.py" execute --manifest "$ADAPTIVE_ROOT/primary225_v1/manifest.json" --reference "$ADAPTIVE_ROOT/original30_pool_v1/REFERENCE.json"
```

This adapter makes native MACE proposals only. It validates exact frozen
reference protocol, optimizer settings, timestamp and exclusion of noncanonical
folds from calibration. Root's `nikasha_pool.py` handles cross-metal native
energies and actual GFN2 solvent terms, then `nikasha_pool_compare.py` performs
the separately calibrated comparison. Missing proposal/solvent terms stay null;
these adapter outputs do not independently contain numerical affinity scores.

Focused real-artifact tests:

```bash
"$ADAPTIVE_PY" -m unittest discover -s tests -p 'test_adaptive_completion*.py' -v
```
