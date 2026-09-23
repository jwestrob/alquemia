# Primary225 transfer submitted

Root authorized full transfer after the independently saved canonical-only
reference became available:25/25 calls, gap2.54669018555 model kcal/mol. That
calibration is development evidence; the forthcoming consumed-fold transfer
remains a separate test of structural robustness.

The initial serial `primary225_v1` is preserved and unlaunched. Current prepared
workspace: `workspaces/adaptive_completion_20260922/primary225_v2_sharded/`.
Manifest SHA256:
`6d94824785d630d294f751b36aa86edbb90ab4cdb7962f8a68658344e5e7ba1c`.

| Fixed task-index shard | Search tasks | Actual submitted job |
|---|---:|---:|
| index modulo4=0 | 103 | 1210021 |
| index modulo4=1 | 103 | 1210022 |
| index modulo4=2 | 102 | 1210023 |
| index modulo4=3 | 102 | 1210024 |

All410 task dictionaries, scientific settings and relevant implementation hashes
exactly match the original unlaunched source manifest. Three focused tests pass
in0.839 seconds, zero skips: exact source replay, one-to-one partition coverage
and invalid-shard rejection/serial compatibility. Full new-manifest dry-run also
passes225/205/410/450 accounting with zero molecular calls. This layout changes
execution only; it adds no scientific candidates, source states or optimizer
starts. The earlier seven numerical/source replay tests remain recorded.

Each job requests1 H200,32 CPUs and200000 MiB. The scheduler determines overlap;
the maximum concurrent request is4 GPUs/128 CPUs/800000 MiB. Actual aggregate
allocation cost is not yet available and will be the sum across jobs, independent
of elapsed overlap. No arbitrary scientific budget or extra runtime cap was added.

Exact submission, selected indices/task IDs, frozen reference and wrapper pins
are in `SUBMISSION.json`. Each job writes a separate
`execution_JOB_shard_N.json` and job-specific GPU log. Collect the complete
225-case/450-endpoint ledger when all four terminate. Preserve all unavailable
and failed results. Root owns the subsequent common-pool scoring and comparisons.
No DFT/GFN2 calls are included in these allocations; production is unchanged.
