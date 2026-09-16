# Derived benchmark ledger

`scripts/ggr_benchmark_export.py` adds the approved development records to the
existing 68-row benchmark. It leaves every inherited field intact, including
historical scores, labels, site vectors, biological groups and decision bands.
The new top-level `development_export` block describes the appended snapshot;
inherited top-level metadata still describes the pinned parent release.

## Executed A/B export

Artifact:
`workspaces/ggr_mechanism_20260915/benchmark_development_ab_v1/benchmark_manifest.json`

- 68 original rows preserved; nine development pairs, 18/18 endpoints complete.
- GGR has six development pairs within its existing `GGR_MglB` group.
- Both alpha-lactalbumin structures stay within one biological observation.
- Aequorin EF3 adds one boundary-control pair; its original six historical
  records remain intact and no site-specific affinity label is added.
- Source-atom/cap graphs, endpoint tasks, unrounded energies, costs, receipts,
  reference records and implementation provenance accompany each addition.
- The inherited alpha preparation records have `protonation_manifest: null`.
  This is reported explicitly; their exact protonated source hashes remain
  pinned. No missing provenance or correction is filled with an invented value.

Command executed from the repository root:

```bash
python scripts/ggr_benchmark_export.py \
  --parent-release workspaces/benchmark_set_20260915/scored_release_1199508/benchmark_manifest.json \
  --collection workspaces/ggr_mechanism_20260915/stage_a_tasks_v1/collection_1199770.json \
  --collection workspaces/ggr_mechanism_20260915/stage_b_tasks_v1/collection_1199802.json \
  --output workspaces/ggr_mechanism_20260915/benchmark_development_ab_v1
```

The existing directory is immutable; a later export requires a new output
directory. Always start from the original parent release, not an earlier
development export. Duplicate stage snapshots, missing rows, altered score
algebra, changed target mappings and inherited classification decisions are
rejected.

## Stage C and conditional C3

The optional repeatable `--sensitivity` argument accepts an actual nominal or
half-step collection. Full collection and task-manifest records, including
gradients, physical maps, endpoints, receipts, bridge checks, failed checks and
unavailable statuses, are retained in GGR's separate
`development_response_records`. The approved GGR 1GLG preparations and plan
must match. A sensitivity collection never creates a displacement benchmark
score, new affinity label, calibrated decision or numerical mechanical
correction. Response and entropy corrections remain null.

A/B completeness and sensitivity completeness are reported separately. The
executed A/B artifact above has no Stage C collection attached; the final
combined artifact should use the final retained C/C3 collections. Failed
physical consistency checks are valid recorded outcomes and do not block
export.

## Verification

```bash
python -m unittest discover -s tests -p test_ggr_benchmark_export.py -v
```

Five tests passed in 6.958 s. These used the actual 68-row release, complete A/B
collections and available executed Stage C artifacts. Corruption tests modify
copies of those records explicitly; they do not supply fabricated successful
scientific outputs. A real-collection copy with one endpoint deliberately
removed verifies unavailable-pair preservation and rejects a false complete
header. The tests launched no calculations. Independent read-only
review found no grouping or protocol-handling blocker. No commit was made by
the exporter author.
