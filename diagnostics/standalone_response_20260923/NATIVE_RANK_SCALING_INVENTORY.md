# Native ORCA GFN2 rank scaling: existing evidence and proposed test

2026-09-23, read-only follow-up requested by root. **No molecular calls, input
changes or submissions were made.** The completed standalone-response branch
remains closed. This question concerns the native ORCA backend used by the
stronger adaptive discriminator, not standalone xTB.

## What is already known

No matched 1/4/8-rank study was found in the inspected project records. A targeted
search of diagnostic plans/reports and operational code found fixed eight-rank
compact workloads and a much larger historical memory failure. Actual receipts:

| Existing workspace | Inspected execution receipts | Recorded ranks |
|---|---:|---:|
| compact_solvation_20260920 | 272 | 8 throughout |
| compact_qualification_20260920 | 48 | 8 throughout |
| slsqp_precision_20260923 | 32 | 8 throughout |
| native_pool_continuation_20260923 | 160 | 8 throughout |

These 512 receipts include different solvers/restart settings and are an inventory,
not a pooled performance benchmark. The 9,141-atom global pilot's 172-rank and
eight-rank attempts produced no converged energy; its memory failure cannot
establish useful scaling for these 100–200-atom contexts. The earlier DDX source
scaling work concerns a different solver.

The recent precision pilot provides a compatible eight-rank cost baseline:
32 native GFN2 cells, individual receipt wall times 13.010724–26.957794s,
median18.4508755s; summed603.464688s (4827.717504 task-rank-seconds). The executor
took89.409918s and its64CPU allocation lasted95s (6080 allocated core-seconds).
These quantities differ from measured CPU utilization; none establishes scaling.
The pinned source cost record retains whole-allocation and actual CPU accounting.

For one real145-atom Q9 Ca/vacuum cell, ORCA printed15.538s total, including
startup1.901s, SCF8.339s and properties2.920s. This suggests both solver and
startup/output costs matter, but is not evidence that any rank count is faster.

## Existing execution constraints

`affordable_workflow.execute` accepts explicit `execution_resources.mpi_ranks`
and `concurrent_tasks`, requiring workers×ranks≤allocated CPUs. The existing
runtime renderer inserts only `%pal nprocs N`; molecular input remains immutable.
`run_orca_task_manifest` runs the actual ORCA binary with `--bind-to none`,
`OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1` and
`OMPI_MCA_hwloc_base_binding_policy=none`, retaining runtime/input/output hashes
and rank receipts. Preserve these fixes and inherited OpenBLAS settings. No need
to change a shared runner or any active225 manifest.

## Proposed smallest matched experiment — not executed

Use the **actual new precision-pilot adaptive_Ca geometry** for Q9Z4J7 (145atoms)
and Q89GY2 (184atoms), each offered identically to Ca and La. Both have Ca charge-3,
La charge-2 and singlet spin. Source tasks and exact coordinates/input/receipt
pins are in [NATIVE_RANK_SCALING_EVIDENCE.json](NATIVE_RANK_SCALING_EVIDENCE.json).

Three fresh rank settings1/4/8 × two sources × two metals × vacuum/ALPB =
**24 calls**. Hold ORCA6.1.1, native mixer, NoAutostart, MaxIter500, electronic300K,
default native parameters, water ALPB, source geometry/state, maxcore2000 and
output settings fixed. Fresh eight-rank cells provide the contemporaneous
reference; preserve old cells separately. Do not seed/restart or choose the
lowest/favorable energy across ranks.

For each configuration run eight workers with CPU-sharing allocations8/32/64CPU
and at least16/64/128GiB respectively, CPU-only on the same established host.
Keep configurations separate in finite manifests and use actual recorded costs;
do not reserve64CPUs for eight serial workers. External load remains a limitation
of a one-batch scaling pilot, so report per-cell latency, summed task-rank-seconds,
batch throughput, allocation cost and memory rather than a universal speed factor.

Predeclare the existing numerical gates: each vacuum/ALPB cell within0.1kcal/mol
of fresh eight-rank, and each source's composite Ca-minus-La contrast within
0.2kcal/mol (MACE is unchanged and reused). Also retain the separate solvent
transfer terms, SCF status/iterations, actual parameter export and runtime rank
confirmation so cancellation cannot hide a failed cell. Any failure or mismatched
state is explicitly unqualified; no automatic retry or backend change. This
numerical study would not recalibrate or alter the running225 science protocol.

Root may authorize this proposed finite comparison separately. This note makes
no claim that reducing ranks is already safe or faster.
