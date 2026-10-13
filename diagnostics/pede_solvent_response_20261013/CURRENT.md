# Active solvent-response investigation — 13 October 2026

Jacob authorized contained work to fix the solvent correction. Read [PLAN.md](PLAN.md).
The saved solvent decomposition completed as3352: all36 origin endpoints close
algebraically. The sample3→1 Ca−La solvent change +21.2566 kcal/mol contains
electrostatics +38.3552, intrinsic remainder −18.6043, hydrogen-bond +1.7026,
surface −0.1969 and reference0. The response is dominated by coupled electronic
and continuum terms, not surface area or a double-counted constant. This does
not prove an unphysical SCF branch or justify deleting any contribution.
Read [the component report](solvent_audit/REPORT.md).

Fixed-origin DFT scout3354 is RUNNING on node-64-768g-5: sample1 La/CPCM water,
native r2SCAN-3c/DefGrid3/TightSCF, twelve MPI ranks. Real preflight passed all
eight frozen input inventories and three archived parser checks. Allocation
48000MiB; the existing runtime derives MaxCore3100MB/rank with25% headroom.
The actual output confirms ORCA6.1.1, native basis and La ECP46/NEL828. No new
energy has completed yet. Only this one molecular cell is submitted. Completion
wakePID349926 and alert-only health watcherPID365231 are armed; commands are in
reference_v1/SCOUT_SUBMISSION.json. Inspect actual result/cost before its vacuum
partner and the remaining six cells. All eight inputs are already pinned.

Cofactor saved audit3353 failed on unavailable-row schema handling; separate
reader recovery3356 completed. All nine admitted origins pass source/PQQ/H/cap
mapping; La1 stays raw-only. All predicted sources have similar long quinone
C–O geometry, including correct controls. Archived crystal states/geometries
are heterogeneous: they do not supply a single justified replacement geometry.
No redox/template substitution or mapped preparation defect was established.
No geometry, score or chemical state changed. Production remains unchanged.

Three saved-artifact jobs used7 allocated CPU-seconds total, including3353's
failure. Molecular scout cost is pending. The first health-monitor attempt
failed on a monitor-only missing resource field; v2 fixes that schema. The
molecular process and original completion wake were unaffected. Both failures
are retained. Do not restart chemistry because a watcher failed.

Exact jobs, receipt paths and wake commands will be recorded in EXECUTION.json.
Before restarting, inspect that record and actual scheduler/output state. Do not
repeat completed work or revive historical LanM/PQQ campaigns.
