> **Final endpoint update, 13 October 2026:** Dy 2666 also converged and passed collection. Both continuations are complete; the six-cell comparison still reverses between Hans sources. Read ../03_comparison_root_recovery/FINAL_REPORT.md. No retry or new molecular work.

> **Update,12 October21:23Pacific:** La2665 completed and its energy/gradients pass the existing collector. Dy2666 remains running; comparison2669 remains queued. Read [the actual La result](LA_RESULT_2665.md). The two-running/pending statements below are historical. No rerun or changed gate.

# Task02 follow-up — Slurm gate passed; two endpoints running

## Executed 12 October 2026, 12:19 Pacific

Real Slurm job **2664 completed 0:0: all eight collector tests passed**, using
one allocated CPU for 29 seconds (pytest 28.11 seconds). This is real archived
output/parser qualification, with zero new molecular calls. Both relocated
positive controls pass with their actual point-charge gradients. Receipt:
`workspaces/claude_lanm_20261012/recovery_02_v1/collector_tests_slurm_v1/RECEIPT.json`.
The earlier six login checks and rejected submissions remain preserved.

The exact two prepared Hans8FNR continuations were submitted only after this
actual pass, with an `afterok:2664` dependency and the receipt enforced again
inside each job. **La job 2665 and Dy job 2666 are RUNNING**, on
`node-64-768g-1` and `node-64-768g-2`, respectively. Both passed the executable
test gate and the unchanged compact force gate. Their runtime sidecars show
12 MPI ranks each, a shared allocation of 48000 MiB each, and MaxCore 3100 MB
per rank derived from that actual grant with the existing 25% headroom.
The current scheduler default is 4000 MiB per CPU; these are not exclusive
whole-node allocations. No memory, walltime, or CPUs-per-task flag was added.

No new seed, method, electronic state, geometry, or convergence gate was changed.
These are the two missing embedded electronic origins, not a new optimization,
a validated affinity calculation, or production promotion. Scientific acceptance
remains pending actual convergence, gradient, seed and embedding checks.

Each job has its own immediate completion wake to root thread
`01a0a63a-ee36-7483-b726-1f7f5b7f75cc`. Watcher PIDs are 2088240 (La) and 2088244
(Dy). Exact submission/wake commands and immutable manifest/test pins are in
each singleton's `SLURM_SUBMISSION_v1.json`. Test watcher PID 2087223 separately
records the real gate completion. The two job receipts are
`La_singleton_origin_v1/completion_event_v1.json` and
`Dy_singleton_origin_v1/completion_event_v1.json` beneath the recovery workspace.
The unchanged `native_scf_health_v5` alert-only watcher also accepts these
singleton manifests directly. It is armed as PID 2089707 (La) and PID 2089708
(Dy), with exact commands in each `HEALTH_WATCH_LAUNCH_v1.json`; it reports
SCF failure/stagnation without cancelling jobs or changing acceptance rules.
Initial real outputs reached point-charge loading and reported that their
first matrix memory checks fit; this is startup evidence only.
The existing EXIT trap collects into `COLLECTION_2665.json` and
`COLLECTION_2666.json`; missing/invalid results remain visible.

Next: inspect those actual terminal collections and Slurm accounting before
recomputing the existing Hans/Mex same-stoichiometry comparison. Do not resubmit
these endpoints merely because the client or watcher stops. Retain both Hans
sources and the four previously accepted DQ2/Mex endpoints.

## Preserved preparation checkpoint, before scheduler recovery


The two fresh Hans8FNR singleton manifests preserve the original scientific
inputs, seeds, electronic states and PBE0-D4 recipe at 12 ranks each. Neither
endpoint has executed. Current scheduler account/association failure blocks
submission; the rejected attempts produced no job IDs.

The original Claude session ran six parser checks on login despite its Slurm
test gate. Those checks passed after preserved failures, but its uncorrupted
relocated controls were invalid because `endpoint.runtime.pcgrad` was absent.
The malformed-copy tests asserted their intended rejection reasons, so this
was an incomplete positive-control harness, not fabricated molecular success.

This follow-up copies the real pcgrad alongside each real fixture and adds both
relocated positive controls. The resulting **eight tests remain UNRUN**. Only
Python AST and shell syntax checks were performed; no imports, preparation,
scientific tests, Slurm submissions or molecular calls were made here.

The singleton wrapper now checks an actual Slurm test receipt before ORCA:
eight exact tests, zero failures/errors/skips, unchanged test/collector/helper
hashes matching the prepared implementation, and test-job `COMPLETED/0:0` from
Slurm accounting. The old login transcript cannot satisfy this gate. It still
requires the unchanged real force gate. No shared executor or prepared artifact
was edited. The current fixed request is 12 tasks, with runtime consistency
checks; it does not derive the scheduler request dynamically.

Both wrappers select healthy `node-64-768g-2`; endpoint submission may explicitly
override to `node-64-768g-1` for the other metal. The runtime guard accepts only
these two reviewed nodes. No memory, CPUs-per-task or walltime flag was added.
The test job requests one CPU for this short serial parser suite; scientific
endpoint jobs retain 12. Qualification of a converged combined MORead/TRAH
restart remains unavailable until the actual endpoint runs.

## Exact next test command after scheduler recovery

From the canonical repository root:

```bash
sbatch \
  "$PWD/workspaces/claude_lanm_20261012/worktree_02/diagnostics/claude_lanm_20261012/02_endpoint_recovery/run_collector_tests.sbatch" \
  "$PWD/workspaces/claude_lanm_20261012/recovery_02_v1/collector_tests_slurm_v1"
```

The new directory must not already exist. Retain the returned test job ID and
arm the established completion watcher. Only after its actual successful gate,
submit each existing singleton with `--dependency=afterok:<actual_test_job_id>`.
Its wrapper now takes two arguments: the unchanged singleton manifest, then
`workspaces/claude_lanm_20261012/recovery_02_v1/collector_tests_slurm_v1/RECEIPT.json`
(both absolute paths). Observe both upstream and endpoint failures and arm
completion wakes. No retry, new seed or method change follows automatically.
