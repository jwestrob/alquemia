# Completed full100 motion-envelope transfer

The primary comparison is complete: **91 correct / 9 unavailable triples**.
Read `REPORT.md` and the saved `COMPARISON_v1.json`; no resubmission is needed.

The frozen `run_v1/INVENTORY.json` records104 unique source/context pairs,
three exact completed pilot pool reuses,101 new pairs and all100 declared triples.
Two source shards contain51 and50 new pairs. Six prior unavailable triples stay in
the denominator. The canonical-only pilot reference is unchanged.

Actual scope:202 native origins,202 bounded searches,at most202 cross-MACE and
404 origin +808 candidate strict nativeGFN2 cells. Native origins completed202/202.
The origin scalar calculations completed normally, but initial path-sensitive collection
rejected the byte-identical runner directory; see `RECOVERY_v2.md`. All original
attempts and zero-search downstream allocations are preserved. Never rerun those
completed origin calculations to repair a collector.

Use the saved origin and downstream `SUBMISSION_*.json` receipts for exact commands
and dependencies. Replacement stage products are under `run_v2`. The private adapter
uses existing MACE/SLSQP/nativeGFN2 kernels; old pinned modules remain unchanged.

```bash
CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
$CPU -m unittest discover -s tests -p test_motion_envelope_transfer.py -v
```

Both terminal pool collections exist. Export an independent comparison without
molecular calls to a **new output name**, preserving the primary result:

```bash
$CPU scripts/motion_envelope_transfer_compare.py \
 --inventory workspaces/motion_envelope_transfer_20260923/run_v1/INVENTORY.json \
 --collections workspaces/motion_envelope_transfer_20260923/run_v2/shard_0/pool/COLLECTION_FINAL.json \
               workspaces/motion_envelope_transfer_20260923/run_v2/shard_1/pool/COLLECTION_FINAL.json \
 --strict workspaces/strict_native_transfer_20260923/run_v1/COMPARISON.json \
 --threefold workspaces/union_triple_transfer_20260923/COMPARISON_v1.json \
 --output workspaces/motion_envelope_transfer_20260923/COMPARISON_replay.json
```

No refit. The exporter retains104 new context/source rows,100 strict triples,
300 declared member occurrences and25 correlated protein groups. Old3.5Å source
contexts are joined per original triple, because multiple old selections can merge
into one envelope selection. Old3.5Å numerical policy remains explicit; the strict
tenfold baseline is the matched scalar comparator. The separate A0A3 Ca3 pilot
abstention is not removed by a favorable La-only transfer result.

## Completed collection recovery

The primary comparison is now complete. The original scheduled collectors failed
on absent `cell_id` metadata. The read-only recovery joins exact case, metal and
candidate fields and verifies physical state; it preserves every original artifact.
Use a **new output name** for an independent replay:

```bash
$CPU scripts/motion_envelope_transfer_collect.py \
 --manifest workspaces/motion_envelope_transfer_20260923/run_v2/shard_0/pool/manifest.json \
 --output workspaces/motion_envelope_transfer_20260923/run_v2/shard_0/pool/COLLECTION_replay.json
```

The saved `COLLECTION_FINAL.json` files and companion `COLLECTION_FINAL_JOIN.json`
are already audited. Prefer them for report-only replays; collection rechecks all
receipts and can take several minutes. No molecular resubmission is needed.
`TESTS_final_v1.txt` records four additional actual-result tests. The final primary
result is 91 correct / 9 unavailable triples. The separate two-seed follow-up owned
by another agent is excluded. `SPREAD.md` records the prespecified matched ranges.
