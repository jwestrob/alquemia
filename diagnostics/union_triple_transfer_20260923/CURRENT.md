# Full100 triple transfer — complete

All molecular tasks and collections are terminal. Final result:92correct/
2inconclusive/6oldunavailable across100 triples, versus94/0/6 for released and
ten-fold precision methods. Both new abstentions concern the Ca-associated A8R3S4
control; no wrong calls.125 source/context pairs give121correct/4inconclusive.
See REPORT.md/RESULT.json and `workspaces/union_triple_transfer_20260923/COMPARISON_v1.json`.
No rollout or refit. Ten actual-fixture tests pass16.704s. Scoped existing-output
A8 audit complete: loss of starting solvent margin dominates; differential
accommodation is slightly smaller, selected motions unchanged for the two failing
members. All four triples and actual3.5026–3.5466Å Ser352-to-Asn300 source contacts
retained. See A8_AUDIT.md. No new molecular tasks in this branch.

The two original CPU stages failed before scalar work due diagnostic-roundoff
comparison. Read-only replay1211145 verified exact selected motions; narrow
validator recovery1211154/55 completed all440 unchanged candidate GFN calls.
Original failures remain in costs:31,153 allocatedCPU-s/644GPU-s. All110 searches
andcrossMACE,38 freshoriginMACE/76originGFN succeeded. Historical running checkpoint
follows for source provenance; do not duplicate its jobs.

Canonical pilot completed and committedb443163:25/25+3/3+separate stress correct.
New independent reference `union_triple_pilot_20260923/REFERENCE_v1.json` is frozen.

Actual source/membership audit:100 predeclared triples,94 prepared/six old missing,
125 unique prepared pairs;70 complete new-ftol pools reused and55 new pools.
All70 actual graph/method/force/receipt audits pass. Largest coordinate-copy
difference1.7763568394002505e−15Å; discrete topology/identities/settings exact.
Failed first read-only audit remains in AUDIT_v1.log. Existing1e−12Å copy tolerance
was applied uniformly to numeric coordinate arrays; no physical bounds changed.

Five actual-fixture prelaunch tests pass11.205s. Origin manifests38MACE/76GFN
validated and jobs1211120/1211121 are complete, all55 paired origins available.
No scientific retries. `run_v2/origins/COLLECTION.json` is the terminal collection.

Two disjoint source shards use the unchanged110-search/cross/440GFN maximum:
shard0 has28sources (GPU1211129, CPU1211130); shard1 has27sources
(GPU1211132, CPU1211133), both validated with zero unavailable inputs.
Check actual submission receipts under
`workspaces/union_triple_transfer_20260923/run_v2/searches/shard_*/` before any
submission. Do not duplicate them or the origin jobs.

Execution helper is immutable `run_v2/execution_v1/union_triple_transfer_run.py`;
dependencies are each task stage's pinned implementation snapshots. Parent owns
production/docs; these are isolated research protocols with noPLM rescore.
Final comparator preserves125 pairs/all100 triples, strict missing members,
both pool policies and released/tenfold contrasts. No reference refit.
