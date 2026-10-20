# Minus endpoint: active CPU work during quiet first SCF iteration

Job15847 completed its initial guess and entered the first incremental Fock/LEAN-SCF
iteration. No completed molecular iteration or convergence residual is available
at this snapshot, and no explicit failure appears.

Batch aggregate CPU time rose05:13:11→05:18:11 over25.03seconds, consistent with
active use of12CPUs. Maximum reported batch RSS35585568KiB stays within48000MiB.
Actual scheduler accounting supports activity; manifest rank estimates were not
used as utilization. Activity alone does not prove useful numerical convergence.

The previously reviewed completed15545 reference took2079.3seconds for its first
iteration. The current quiet phase is compatible with that same-recipe history;
its independent source review is linked through ROOT_HEALTH_REVIEW.json.

**Continue unchanged.** No restart, cancellation, duplicate calculation, changed
SCF/grid/guess or watcher modification. Native and physical admission remain pending.
