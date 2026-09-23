# Completed full225 operations and read-only replay

All13 finite allocations/collectors completed. Do not resubmit molecular jobs
or alter their optimizer settings. The saved COMPARISON_v1.json is the final
225-source result; TRANSFER225_REPORT.md and TRANSFER225_RESULT.json summarize it.

## Fixed chain

| Shard | Sources / searches | Proposal | Prepare + GFN2 | Cross-MACE + collection |
|---|---|---:|---:|---:|
|0|51 /102|1210573|1210586|1210587|
|1|51 /102|1210574|1210588|1210589|
|2|51 /102|1210575|1210590|1210591|
|3|51 /102|1210576|1210592|1210593|

All edges use `afterany`, retaining actual failed endpoints/cells. The existing
CPU/GPU executors only run valid prepared tasks. The final comparison job1210601
depends on all four collection jobs. Exact commands/wrappers/pins are under
`workspaces/union_adaptive_20260923/transfer225_v1/`; production is unchanged.
The four solver jobs were subsequently broadened in place to
`gpu,standard-shared`;64CPUs/128GiB, job IDs, priorities and dependencies stayed
unchanged. `CPU_PARTITION_UPDATE_v1.json` preserves each before/after scheduler
record. The original wrapper's submission directive remains historical evidence.

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
UNION_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
UNION_STAGE=$PWD/workspaces/union_adaptive_20260923/transfer225_v1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sacct -X -j 1210573,1210574,1210575,1210576,1210586,1210587,1210588,1210589,1210590,1210591,1210592,1210593,1210601 -o JobID,State,Elapsed,AllocCPUS
"$UNION_PY" scripts/union_adaptive_transfer.py validate_selection \
  --inputs "$UNION_STAGE/INPUTS_shard_0.json" \
  --calibration workspaces/consistent_context_20260922/calibration28_v1/collection_final_v2.json
"$UNION_PY" -m unittest discover -s tests -p test_union_adaptive_transfer.py -v
```

All four terminal collections exist. The comparison can be replayed with
no molecular work into a new output path:

```bash
"$UNION_PY" scripts/union_adaptive_transfer_compare.py \
  --selection "$UNION_STAGE/SELECTION.json" \
  --collections "$UNION_STAGE/shard_0/pool/collection_final.json" \
                "$UNION_STAGE/shard_1/pool/collection_final.json" \
                "$UNION_STAGE/shard_2/pool/collection_final.json" \
                "$UNION_STAGE/shard_3/pool/collection_final.json" \
  --baseline workspaces/nikasha_recovery_20260922/proposal_comparison.json \
  --union_static workspaces/consistent_context_20260922/primary225_v1/comparison_v1.json \
  --native workspaces/adaptive_minimal_pool_20260923/COMPARISON_recovered_v1.json \
  --standalone workspaces/standalone_xtb_transfer_20260923/run_v1/COMPARISON_v1.json \
  --output "$UNION_STAGE/COMPARISON_replay.json"
```

All225 sources,17 existing exclusions, four reused pilot folds and every new
search/solver failure remain visible. Groups use strict La4/Ca5 medians and equal
balanced means; all100 declared La triples stay in the denominator. Each method
uses its already frozen reference; this operation cannot refit thresholds or
select replacement members. New scalar values are not binding free energies.

Input construction is the recorded `union_adaptive_transfer.py selection`
operation with explicit transfer/sources/calibration/reference/pilot/agreement
paths in `SELECTION.json`. Each `INPUTS_shard_i.json` then uses the existing
`union_adaptive.py prepare` interface. The molecular preparation commands and
measured `/usr/bin/time` receipts are preserved as `PREPARATION_shard_i.log`.
