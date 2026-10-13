# Active solvent-response investigation — 13 October 2026

Jacob authorized contained work to fix the solvent correction. Read [PLAN.md](PLAN.md).
The saved solvent decomposition completed as3352: all36 origin endpoints close
algebraically. The sample3→1 Ca−La solvent change +21.2566 kcal/mol contains
electrostatics +38.3552, intrinsic remainder −18.6043, hydrogen-bond +1.7026,
surface −0.1969 and reference0. The response is dominated by coupled electronic
and continuum terms, not surface area or a double-counted constant. This does
not prove an unphysical SCF branch or justify deleting any contribution.
Read [the component report](solvent_audit/REPORT.md).

Fixed-origin DFT scout3354 converged: sample1 La/CPCM water, native
r2SCAN-3c/DefGrid3/TightSCF. Slurm reports FAILED because the collector looked
for the final TRAH tolerance in an optional absent summary row. Actual printed
threshold1e-5 and final macro residual5.663e-7 pass the unchanged criterion.
Collection-only recovery3370 passed15 real-artifact/parser checks with zero
molecular reruns. One of eight endpoints is now available. Read
[SCOUT_REPORT.md](SCOUT_REPORT.md) and the versioned recovery report. The
original failed batch/collection remains preserved; no tolerance was loosened.

The remaining seven frozen endpoints are RUNNING in four disjoint12-MPI lanes
on node-64-768g-5:3371 sample1 La vacuum;3372 sample1 Ca water/vacuum;3373 sample3
La water/vacuum;3374 sample3 Ca water/vacuum. All four started without a launch
error. `REFERENCE_DISPATCH_v1.json` maps the eight tasks without repeating the
scout. Each pair stops if its first cell is unavailable. No optimization or new
chemical state. Completion and alert-only health wakes are armed for every lane;
exact commands/PIDs are in the linked CONTINUATION_SUBMISSION.json files and
reference_continuation_v2/SUBMISSION.json. Shared QoS respected,48 new CPUs total.

Actual scout cost9580s on12 CPUs =114960 allocated CPU-seconds; sampled
peakRSS18172MiB. Actual memory grant48000MiB, runtime MaxCore3100MB/rank with25%
headroom. This is development reference cost; no routine-method speed claim.
No solvent transfer or source double difference exists until matching cells
complete. Root owns execution; comparison agent prepares the saved-result join.

Cofactor saved audit3353 failed on unavailable-row schema handling; separate
reader recovery3356 completed. All nine admitted origins pass source/PQQ/H/cap
mapping; La1 stays raw-only. All predicted sources have similar long quinone
C–O geometry, including correct controls. Archived crystal states/geometries
are heterogeneous: they do not supply a single justified replacement geometry.
No redox/template substitution or mapped preparation defect was established.
No geometry, score or chemical state changed. Production remains unchanged.

Three saved-artifact jobs used7 allocated CPU-seconds total, including3353's
failure; parser recovery3370 adds2. The first health-monitor attempt
failed on a monitor-only missing resource field; v2 fixes that schema. The
molecular process and original completion wake were unaffected. Both failures
are retained. Do not restart chemistry because a watcher failed.

Exact jobs, receipt paths and wake commands will be recorded in EXECUTION.json.
Before restarting, inspect that record and actual scheduler/output state. Do not
repeat completed work or revive historical LanM/PQQ campaigns.
