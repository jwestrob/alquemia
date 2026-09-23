# Full100 motion-envelope transfer

The frozen `run_v1/INVENTORY.json` records104 unique source/context pairs,
three exact completed pilot pool reuses,101 new pairs and all100 declared triples.
Two source shards contain51 and50 new pairs. Six prior unavailable triples stay in
the denominator. The canonical-only pilot reference is unchanged.

Actual scope:202 native origins,202 bounded searches,at most202 cross-MACE and
404 origin +808 candidate strict nativeGFN2 cells. Native origins completed202/202.
The scalar calculations completed normally, but initial path-sensitive collection
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

Once both terminal pool collections exist, export a new comparison without calls:

```bash
$CPU scripts/motion_envelope_transfer_compare.py \
 --inventory workspaces/motion_envelope_transfer_20260923/run_v1/INVENTORY.json \
 --collections workspaces/motion_envelope_transfer_20260923/run_v2/shard_0/pool/COLLECTION_FINAL.json \
               workspaces/motion_envelope_transfer_20260923/run_v2/shard_1/pool/COLLECTION_FINAL.json \
 --strict workspaces/strict_native_transfer_20260923/run_v1/COMPARISON.json \
 --threefold workspaces/union_triple_transfer_20260923/COMPARISON_v1.json \
 --output workspaces/motion_envelope_transfer_20260923/COMPARISON_v1.json
```

No refit. The exporter retains104 new context/source rows,100 strict triples,
300 declared member occurrences and25 correlated protein groups. Old3.5Å source
contexts are joined per original triple, because multiple old selections can merge
into one envelope selection. Old3.5Å numerical policy remains explicit; the strict
tenfold baseline is the matched scalar comparator. The separate A0A3 Ca3 pilot
abstention is not removed by a favorable La-only transfer result.
