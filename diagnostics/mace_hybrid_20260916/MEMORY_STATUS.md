# MACE memory-work status

**Completed 2026-09-16.** Job 1200381 completed all 12 calls on one RTX A5000.
Both 9,141-atom primary endpoints took about 58 seconds each, with 9.55 GiB
peak GPU allocation. Maximum host RSS across all full calls was 1.65 GiB.
No host offload was needed. No memory-work job or H200 continuation remains active.

All 17 tests and the four original-core comparisons pass. Full-system repeats,
translations and charge closure pass. Rotation and the hybrid partition check
fail their frozen tolerances. Baseline/default, weights and physical inputs
remain unchanged; large-model compatibility is untested.

See [final results and exact costs](MEMORY_RESULTS.md), [runbook](RUNBOOK.md),
[agreement](MEMORY_AGREEMENT.md) and [implementation](BLOCKED_KERNEL.md).
Current campaign: `workspaces/mace_hybrid_20260916/blocked_v6`.
Terminal collection: `collection_job_1200381.json`; report: `REPORT.md`.
The separate accounting watcher writes `job_1200381_accounting.json`.

Preserved memory-development allocations: 1200365, 1200368, 1200369, 1200370,
1200371, 1200372, 1200379, 1200380 and 1200381. Own pending H200 job 1200309
was cancelled without compute after initial core verification. Offload job
1200372 had already failed when cancellation was requested after excessive RSS
was observed; the correction receipt records this precisely.
