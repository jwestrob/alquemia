# Reporting overhead: isolated engineering check

The unchanged timing campaign remains authoritative for end-to-end speed.
No model, geometry, label, threshold, endpoint or production default changes here.

Profiling the already completed first canonical case found1,948XYZ parses in a
single report. Calibration verification alone accounted for220.17of285.69profiled
seconds; the whole replay consumed276.91CPU seconds. Profiling overhead is large:
these instrumented times must not be compared directly with ordinary timings.
The replay's complete result file exactly matched the original report hash.

Isolated `interface_source_v5` extends the existing operation-local file cache
to parsed XYZ coordinates. The strict parser and float conversion are unchanged.
Every lookup still checks file identity/mtime/ctime with NFS revalidation; recent
files are not reused. Cached atoms are immutable tuples and each caller receives
a fresh list. Production source and the runningV2 snapshot are untouched.
The exact change fromV4 is retained in XYZ_CACHE.patch.

Four real-artifact tests pass in18.552s: all54supported PQQ endpoint coordinates
are bitwise equal; parsing is reused only within one operation; caller mutations
do not change saved coordinates; corrupted/deleted real fixture copies are not
served from cache. No invented energies or successful scientific calculations.

An unprofiled local cached report passed with identical scientific fields;
only `report_implementation` differed. It took101.3513wall/93.2581CPU seconds on
the login host. This cannot establish speed versus the65.49s GPU-host report.

For a matched engineering measurement, job1201589 executes exactly two report
invocations on the original CPU host, V4uncached thenV5cached, on the same already
computed first-case inputs. One CPU,4096MiB,no GPU,new model evaluations=0.
Full commands/source pins are in
`workspaces/mace_pqq_utility_20260918/report_cache_matched_v1/manifest.json`.
Acceptance: identical scientific fields and lower measured report wall/CPU time.
No combined end-to-end speed claim may be inferred from these stage measurements.

## Matched result: passed

Job1201589 completed. Both reports have identical scientific fields; only the
recorded reporter implementation path differs. Uncached reporting took62.06s;
cached reporting took44.67s (**1.3893x**,17.39s saved). CPU time fell from53.71to
36.80s. Maximum process RSS was584896versus593828KiB. Allocation cost107core-s,
zero GPU time and zero model evaluations. Exact GNU-time outputs, Slurm accounting,
report hashes and joined result are retained in report_cache_matched_v1.

This is a measured reporting-stage improvement. The final end-to-end comparison
will use the unchanged full25case timing campaign; no unmeasured combined latency
or additional successful model calls are implied. V5 remains isolated and opt-in.
