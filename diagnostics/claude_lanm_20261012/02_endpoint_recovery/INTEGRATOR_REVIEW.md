# Task02 follow-up — prepared, tests UNRUN

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
