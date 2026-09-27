# Hans EF3 scout: Dy nonconvergence stops this attempt

The four-endpoint attempt is incomplete. Both La endpoints completed; both Dy
calculations were cancelled after persistent severe SCF instabilities. No La/Dy
response contrast, affinity or discrimination result is available. Do not run the
two prepared transfer sources or retry this unchanged setup.

## Retained results

| Endpoint | Status | Energy, Hartree | Torsion gradient, kcal/mol/radian |
|---|---|---:|---:|
| La/A | Complete | -5730.072520781592 | -9.525688543235155 |
| La/B | Complete | -5730.072802298282 | -1.5827529285371433 |
| Dy/A | Attempted, unconverged at cancellation | unavailable | unavailable |
| Dy/B | Attempted, unconverged at cancellation | unavailable | unavailable |

La work for the declared +2degree Asp85 motion is -0.17665439002530572 kcal/mol.
That is a local electronic response only. It does not establish metal specificity,
an equilibrium structure, complete protein mechanics or qualified hybrid forces.
La/A and La/B took approximately2h48m and2h6m. The archived partial Dy logs show
large TRAH residuals and negative orbital-gap warnings, not a nearly converged
endpoint. These observations do not uniquely establish the underlying cause.

The generic collector labels Dy rows `missing` because cancellation prevented
completed execution receipts. Preserve that immutable collection, but distinguish
this from unattempted chemistry: both calculations ran and their partial outputs
are retained. EARLY_STOP_20260927.json preserves evidence at cancellation.

## Actual cost and monitoring failure

Job1219501: CANCELLED, 29,872seconds ×224 allocatedCPUs =6,691,328 allocatedCPU-s
(1,858.7022 core-hours). Collector1219502:3seconds ×1CPU. Total6,691,331 allocated
CPU-s, zeroGPUs. Batch-step MaxRSS168,090,916KiB is the scheduler's reported step
memory statistic, not a claimed measured node-wide peak or CPU utilization.

After both La tasks ended, only two56-rank Dy tasks remained, while the allocation
continued holding224cores. The completion watcher successfully delivered the
terminal event, but it did not monitor SCF progress or idle allocation inside the
job. This monitoring omission allowed prolonged unproductive execution. A runtime
convergence/stagnation monitor and a smaller first electronic-state qualification
are prerequisites to any further expensive attempt; neither has been implemented
or launched by this collection step.

## Stop and recovery state

Root cancelled only owned worker1219501 and left its afterany collector intact.
Both are terminal; the completion event is acknowledged in ROOT_ACK.json. No
other session's jobs were touched. All source preparations, outputs and completed
La receipts remain available under lanm_ef3_hans_scout_v1. Prepared Hans8FNR/Mex8FNS
remain unscored. No production change, automatic retry, new submission or push.

Next action is scientific diagnosis of the saved Dy histories and model choice,
not increasing iteration limits or repeating this calculation unchanged.
