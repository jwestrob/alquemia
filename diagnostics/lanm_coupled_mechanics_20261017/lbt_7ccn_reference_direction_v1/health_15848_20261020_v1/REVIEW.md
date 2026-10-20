# Quiet first SCF iteration; continue unchanged

Job15848 is running its first incremental Fock/LEAN-SCF iteration. Initial guess
completed; no completed molecular iteration or convergence residual is available
yet. The native snapshot contains no explicit failure.

Two scheduler samples show aggregate batch CPU time increasing from04:51:15 to
04:56:15 over25.03seconds, consistent with active use of12CPUs. Maximum observed
RSS is35582868KiB within the48000MiB allocation. This is actual accounting;
manifest rank estimates were not treated as utilization.

The completed same-recipe15545 reference printed2079.3seconds for its first DIIS
iteration and2082.7seconds for a later incremental Fock rebuild. Arendt independently
reviewed these source lines (2354–2359 and2368–2369). The current silent phase is
compatible with long native computation. CPU activity alone does not establish
useful convergence, and no new numerical or physical admission is claimed.

**Decision: continue unchanged.** No restart, cancellation, additional calculation,
SCF/grid/guess change or watcher modification. Original finite-test ownership,
completion wakes and alert-only health policy remain active.
