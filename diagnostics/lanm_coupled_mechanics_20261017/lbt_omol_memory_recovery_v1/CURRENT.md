# Resource-only central recovery — running

21 October 2026 UTC. **Job 16680 is running.** Actual startup confirms 12 MPI
processes, one numerical thread each, inside 24 allocated CPUs and the same
96,000 MiB default-memory entitlement. MaxCore is 6,200 MB per process.
The molecular template, XYZ and basis are byte-identical to failed job 16611;
the runtime input changes only PAL24→12 and MaxCore3100→6200.

Preparation 16667 passed in 7 CPU-seconds with zero molecular calls. Result
cd2d06b2 / manifest d0c662ef / archive db5e34b1 preserve the source and input
checks. Independent actual review a0e929b9 and root review bee284c2 passed.
Source 5cad2623, specification 29b32616, run configuration 908678ef and root
execution decision c996c2e2 bind this one endpoint.

Root startup review 5fdad85d confirms the grant/runtime and native SCNL,
libXC and La ECP46 headers. Completion watcher 91376 and alert-only health
watcher 91377 are active. PModel initialization and beta rotation have now completed. The final molecular
energy/gradient are pending; no force or state qualification is established.
Independent startup review eae72570 agrees. Run submission archive 6ec11a5c
preserves immutable inputs, configuration, release, runtime and watcher records.
It is a startup snapshot, not a completed molecular result.

## Interim health review,04:31UTC

The silence alert was reviewed: all12 LEANSCF workers accrue aboutone core of
CPU time each; cgroup memory is67.7GiB against93.75GiB, zero swap/OOM events.
Output remains at the first Fock build with no SCF iteration or residual yet.
Continue16680 unchanged; no duplicate or cancellation. See
health_review_16680_silent_v1/REPORT.md and its exact snapshots.

## Why this recovery

The prior 16611 PModel initialization succeeded, but its first molecular Fock
calculation exceeded memory: ORCA estimated 4,650.8 MB per process against
3,100 MB available. It produced no final energy or gradient. Its failure,
23,088 CPU-seconds and native evidence are preserved in the parent PModel
branch, archive 94344987. Original 16602 guess failure is also preserved.

This finite attempt changes only the resource layout. No memory/time override,
geometry, Hamiltonian, initial guess, state, threshold or admission-rule change.
Fewer ranks may change numerical reductions and the SCF trajectory; sufficient
memory and convergence remain unproven. Root owns integration. No duplicate,
automatic retry, displaced pair or downstream molecular calculation is released.
