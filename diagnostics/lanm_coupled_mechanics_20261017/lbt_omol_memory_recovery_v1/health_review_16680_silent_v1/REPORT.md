# Job 16680: silence alert reviewed

21 October 2026,04:30–04:31UTC. **Continue unchanged.**

Two host-process observations10.009seconds apart show all12 owned LEANSCF
workers running and each accruing approximately10.00–10.02CPU-seconds.
The main output is quiet during the first incremental Fock build, after
successful PModel initialization and beta rotation. No first SCF iteration,
energy/residual row or gradient is available yet.

The job cgroup uses72,659,742,720bytes against100,663,296,000bytes allowed
(about67.7/93.75GiB). Peak72,887,996,416bytes; zero swap, OOM or OOM-kill events.
The previous3100MB-versus4650.8MB warning is absent; native stderr is empty.
The extra12 allocated CPUs are the declared memory entitlement, not12 missing
MPI workers. Slurm accounting and cgroup memory measure different quantities;
process RSS values must not be summed as a unique physical-memory measure.

The alert establishes output silence, not a stalled calculation. CPU activity
and current memory headroom support continuation, but do not establish useful
SCF convergence or future memory sufficiency. No cancellation, duplicate run,
parameter change or new scientific calculation follows. Existing completion
and alert-only health notifications remain active.

ROOT_REVIEW.json binds raw process, cgroup, scheduler, receipt and output
snapshots. Live watcher files are unchanged. This is an interim observation,
not a terminal result or final resource-cost report.
