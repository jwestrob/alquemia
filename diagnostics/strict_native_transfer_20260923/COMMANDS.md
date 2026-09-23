# Exact finite commands

Run from the repository root with the existing pinned Python:

```bash
PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
RUN=workspaces/strict_native_transfer_20260923/run_v1
```

Preparation and all four independent dry-runs completed in `PREFLIGHT.json`.
All jobs1211337–1211340 use the corresponding snapshot in `RUN/implementation`;
`SUBMISSION.json` records the exact commands and resource allocation. Do not
resubmit these existing tasks. No restart seeds are required for fresh scoring.

Inspect or collect one existing shard (change only the explicit shard index):

```bash
$PY "$RUN/implementation/strict_native_transfer.py" validate \
  --manifest "$RUN/shard_0/manifest.json"
$PY "$RUN/implementation/strict_native_transfer.py" collect \
  --manifest "$RUN/shard_0/manifest.json" \
  --output "$RUN/shard_0/COLLECTION_recovered.json"
```

The submitted wrapper collects `COLLECTION.json` even when execution fails and
preserves the executor exit code. Existing failed attempts remain visible;
there is no automatic scientific retry. Use a new collection filename for a
report-only recovery; never overwrite old records or rerun the molecules.

Once all four actual collections exist, compare the full panel without new
molecular calls or calibration:

```bash
$PY "$RUN/implementation/strict_native_transfer.py" compare \
  --inventory "$RUN/INVENTORY.json" \
  --collections "$RUN/shard_0/COLLECTION.json" "$RUN/shard_1/COLLECTION.json" \
                "$RUN/shard_2/COLLECTION.json" "$RUN/shard_3/COLLECTION.json" \
  --output "$RUN/COMPARISON.json"
```

Eight reused pools are loaded from the pinned strict32 fresh collection and
water_basins' separately owned four-pool comparator. Every missing cell remains
unavailable; the old precision/DFT/static fields remain unchanged beside the
new strict results. Strict group medians require every declared member.

Real-fixture tests:

```bash
$PY -m unittest discover -s tests -p test_strict_native_transfer.py -v
```

Fresh preparation for a separately named reproduction directory requires the
explicit inputs in `PLAN.md` and the already completed comparator manifest;
it is not needed to recover the current jobs/results.
