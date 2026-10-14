# ExaA donor-response checkpoint

**Two native DFT scouts are running:4109 origin and4110 adaptive_La.** Each uses
12MPI ranks on the shared standard node, normal priority, with completion and
health wakes armed. Collector4111 waits afterany on both and keeps all eight
planned cells visible, including six not yet submitted. No completed new energy
is claimed at this checkpoint. Native output confirms the expected La46-core ECP.

Jacob authorized autonomous next rounds on14October2026. Read PLAN.md. This test
asks whether independent vacuum DFT supports the large MACE stabilization on a
saved ExaA donor path. It preserves the exact192-atom origin/adaptive_La geometry,
charge, cofactor state and caps. The matched control is sample0; target is sample3.
This common-path test is distinct from full operational candidate minimization.

## Verified execution

- Slurm4108 completed in6s/12CPUs: actual geometry/state/electron/mapping gates
  and all11 real archived parser/algebra checks passed; zero molecular calls.
- Ca has838 explicit electrons; La has874 total/828 explicit after ECP46. Both
  are the declared singlets. Existing dry oxidized-PQQ3− hypothesis retained.
- Actual shared grant48,000MiB; renderer derives MaxCore3100decimalMB per rank
  with existing headroom. No memory/time override or scheduler priority boost.
- Initial native endpoints4109/4110 are active, each at the expected ECP stage
  when inspected. No SCF convergence or response result yet.
- Read LAUNCH.json and immutable workspace SUBMISSION.json files for exact pins,
  commands and watcher PIDs. Only root owns submissions/collection/integration.

## Resume precisely

Workspace: `workspaces/exaa_donor_response_20261014/reference_v1/`. Manifest SHA:
`d37733923c3203f41b73d5a3c1789d396e22f8da69e624df2c2d724d0141bcf1`.

Inspect jobs4109/4110/4111 and task receipts before any action. Each native job
writes `attempts/JOB_execute_TASK/COLLECTION.json`. Joined actual scout results:
`workspaces/exaa_donor_response_20261014/scout_comparison_v1/execution/RESULT.json`.
Absent outputs remain absent. Re-arm recorded watchers after host restart; never
resubmit chemistry simply because a watcher died. Delayed4108 preflight event was
already read and acknowledged while root was active; no separate idle-wake claim.

Review both scout outputs for admitted energy, final solver residual, orbital
warnings/state ambiguity and measured cost. The first meaningful quantity is
La E(adaptive_La)−E(origin). Then pursue the remaining six declared endpoints if
this reference is usable. All source/implementation agents completed their work.
No extra approval is required within the current authorized scope.

All three current job IDs and actual watcher commands are in LAUNCH.json-linked
receipts. The eight-cell matrix is not complete; production/calibrations and old
negative results remain untouched. No new folds, optimization or parameter fit.
