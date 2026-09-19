# Replay and prepare the scanner comparison

Run from the repository root. Outputs are write-once; choose a fresh replay path.
Job1202084 is complete; do not resubmit it.

```bash
SCANNER_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$SCANNER_PY" scripts/hydration_scanner.py collect \
  --manifest workspaces/hydration_scanner_20260919/mace_v2/manifest.json \
  --output workspaces/hydration_scanner_20260919/recollection_v1.json
"$SCANNER_PY" -m unittest discover -s tests -p test_hydration_scanner.py -v
```

Independent preparation, no new energies:

```bash
"$SCANNER_PY" scripts/hydration_scanner.py prepare \
  --inventory workspaces/site_classifier_20260918/inventory_v2/inventory.json \
  --core-collection workspaces/hydration_network_20260918/core_transfer_v1/collection_1201867.json \
  --normalized-collection workspaces/hydration_square_20260918/repaired_v1/collection_1201801.json \
  --ggr-comparison diagnostics/hydration_network_20260918/GGR_REPLICATE_RESULT.json \
  --agreement diagnostics/parallel_discriminator_20260919/SCANNER_COMPARISON.md \
  --output workspaces/hydration_scanner_20260919/reprepare_v1
"$SCANNER_PY" workspaces/hydration_scanner_20260919/reprepare_v1/implementation/mace_hybrid.py \
  dry-run --manifest workspaces/hydration_scanner_20260919/reprepare_v1/manifest.json
```

The actual allocated execution command is preserved in
`workspaces/hydration_scanner_20260919/mace_v2/submission.json` with its necessary
QOS-duration correction in`scheduler_duration_fix.json`. It uses the established
`diagnostics/mace_hybrid_20260916/run_pilot.sbatch` and frozenworkerimplementation.
No new workflowdaemon, fallbackscorer, or altered scientificdefault is installed.
