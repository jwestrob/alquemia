# V3 saved component-replay diagnostic

This is a diagnostic version of unchanged ANALYZE_V2, following saved-only job13455.
All four original cells failed exact raw-component/full-gradient replay. No cause
is assigned from the failure alone. Preserve both prior failures and all gates.

Exactly the same four frozen POINTs, reference, units, archived Jacobians and
component order are read. Before the original exact replay require, record each
component's exact-equality flag, maximum absolute/signed difference and its
physical atom/axis, number and indices/IDs of differing entries, and stored versus
replayed values. Record stored-component sum versus stored full gradient,
raw-replayed sum versus stored full gradient, and raw versus stored component sums.
These diagnostics apply no tolerance and do not turn a failure into admission.

Record Python/NumPy, build/BLAS configuration, hostname/architecture, thread and
Slurm CPU environment. Preserve component/map pins and cap-anchor identities.
Do not recompute a map, change scientific settings or claim host/BLAS roundoff
without the actual measurements. Any remaining mismatch retains the original
unavailable result and exit failure; diagnostic files remain linked in RESULT.

One CPU under Slurm, original NumPy2.4.6 environment and one thread. No native
process, model load, molecular call, geometry, parameter change or login numerical
work. No new test harness or retry is introduced. Root runs this once after review.
