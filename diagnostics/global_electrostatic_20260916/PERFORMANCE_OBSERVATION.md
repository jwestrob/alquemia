# TABI performance observation — 2026-09-16

Read-only inspection at **01:14–01:18 PDT**, while jobs **1199956** and
**1199959** ran on `node-48-256g-13`. No jobs, benchmarks, profiler attachments,
affinity changes, solver changes, or reruns were performed for this inspection.
These are execution observations, not completed scientific results or a runtime
forecast.

## Verified resource layout

`lscpu -e`, `/proc/PID/status`, and `ps` establish:

| Item | Observation |
| --- | --- |
| Hardware | Two Intel Xeon E5-2695 v2 sockets; 12 physical cores/socket; 2 hardware threads/core; 48 logical CPUs total |
| Initial pair, job 1199956 | Both processes allowed CPUs `0,24`; these are **siblings of physical core 0**, socket 0 |
| Later job 1199959 | All processes allowed odd logical CPUs `1,3,...,47`; these cover **12 physical cores on socket 1** |
| Active work | Initial 2 plus later 21 native solves; 2 later isolated solves had already finished |
| Native threading | Every active `tabipb` process had exactly 1 thread |
| CPU use | Initial processes about 99.6%; later processes about 94–96%; nice values all 0 |
| Memory | About 252,855 MiB available on the node; native RSS about 65–138 MiB/process |
| Pressure | Current CPU, memory, and I/O pressure averages were 0.00 at inspection |

Thus the reported logical CPU allocation is not a count of independent physical
cores. The first pair shares one physical core; the later 21 solves share 12
cores. High `%CPU` reports scheduled hardware-thread time and does not imply an
unshared core or a particular clock rate.

## Clock observations and their limits

At 01:15:26, `lscpu -e=CPU,CORE,SOCKET,MHZ` reported busy CPUs 0/24 at
2433.912/2437.301 MHz. Busy socket-1 CPUs mostly reported approximately
872–906 MHz. A separate read of `scaling_cur_freq` at 01:15:49 returned:

| Logical CPU | Reported kHz |
| --- | ---: |
| 0 | 2400010 |
| 24 | 2400010 |
| 1 | 909690 |
| 25 | 891321 |

All four used `schedutil`, with configured min/max 1200000/3200000 kHz. The
below-minimum readings deserve operational investigation, but **their cause is
not established here**. Selected core/package thermal-throttle counts were
zero; package temperatures were 41/38 °C, critical temperatures 87 °C, and
critical alarms zero. These observations do not support a thermal-throttling
claim. No frequency/governor setting was changed.

A raw two-snapshot capture is preserved at
`workspaces/global_electrostatic_20260916/performance_observation/node-48-256g-13_20260916.json`,
SHA256 `d2b32e30cb8881b131296451093349c8af12bd3ce4c9e682ff6c4ac5e87fd100`.
It includes `lscpu`, topology, `/proc/cpuinfo`, `/proc/stat`, per-process
`stat/status/schedstat/cwd`, all readable CPU frequency and thermal-throttle
files, hardware-monitor temperatures, memory, and pressure counters. Snapshots
at 01:17:23.037402 and 01:17:33.659648 were 10.62223856151104 s apart, with
10 s of passive waiting between reads; no workload was launched.

Every active solver accumulated approximately **10.60–10.63 CPU seconds** in
that interval. Socket-1 frequency medians across its logical CPUs were
923.092 and 907.1875 MHz, while CPUs 0/24 read 2110.646–2663.916 MHz across
the snapshots. All captured thermal-throttle counters remained zero; package
temperatures were 42 °C and 37–38 °C. This confirms ongoing execution and a
persistent reported frequency disparity. It does not establish its mechanism.

Representative `/proc/PID/schedstat` readings showed negligible run-queue delay
relative to runtime (initial La approximately 1.35 s versus 1708.84 s; repeated
La approximately 0.015 s versus 1067.89 s). Representative NUMA maps placed
large anonymous arrays on the process's local socket. This does not measure
cache or memory-bandwidth contention.

## Same-work slowdown and different-work checks

Existing immutable NanoShaper receipts establish a same-mesh slowdown:

- Initial `1h4i_qm33_La/primary`: **4.01167456060648 s**.
- Later `1h4i_qm33_La/repeat`: **34.7172104306519 s**.
- Both have 109,054 vertices and 218,056 faces and the same prescribed mesh.

Receipt paths are
`workspaces/global_electrostatic_20260916/surfaces_v1/1h4i_qm33_La/{primary,repeat}/attempt_0001/mesh/execution.json`.
This is mesher elapsed time under different concurrent load, not a controlled
benchmark or an established native-solver slowdown factor.

The numerical variants also differ in intrinsic work:

- Refined meshes have **246,760 vertices / 493,468 faces**, about 2.26 times the
  primary surface size; their matrix-vector operations are not directly
  comparable with primary ones.
- The tree check keeps the primary mesh but raises degree 5 to 7 and lowers
  opening parameter 0.5 to 0.3. The pinned native source contains three nested
  interpolation loops for particle/cluster interactions and six for
  cluster/cluster interactions; higher degree changes the cost substantially.
- `particles.cpp:187` computes the source term with surface-by-charge loops.
  `boundary_element.cpp:60` performs the repeated matrix-vector work.
  `boundary_element.cpp:655` then computes energies after GMRES finishes.
  Detailed native timings become available only at completion.

Pinned source:
`workspaces/global_electrostatic_20260916/software/source_inspection/`,
TABI commit `fe1c237b057418fed48535db125394607040d9de`.
The initial pair had reported 9 GMRES iterations near 27 minutes, while later
primary-equivalent tasks had generally reported 1–2 near 17 minutes. Those
checkpoints are not equal-work elapsed-time measurements and give no justified
completion forecast.

## Practical conclusion

There is verified hardware-thread sharing and a large contemporaneous socket
clock disparity, plus more expensive refined/tree tasks. The full contribution
of each effect is unresolved; cache contention, power management, and other
causes must not be asserted from these observations alone. Memory exhaustion,
I/O stalls, or an idle solver are not supported by the inspected state.

For future runs, request and verify one physical core per native solve and
record actual affinity/clock behavior before using allocated logical CPU-hours
as a hardware-independent affordability measure. This is a recommendation;
the current jobs and frozen numerical schedule were left untouched. Final
per-task timers and execution receipts remain the authoritative measured cost.

If the parent elects technical recovery of the affected tasks on a verified
healthy allocation, retain their frozen inputs, numerical settings, existing
attempts, and all consumed cost. The initial pair on socket 0 need not be
disturbed. No recovery was executed by this inspection.

## Recovery observation — 01:20–01:21 PDT

The parent subsequently canceled only affected job 1199959 and submitted the
same frozen tasks as array 1199964, excluding node-48-256g-13 and requesting
`--hint=nomultithread`. The initial pair continued. This subsection inspects
the recovery; it did not submit or modify it.

The first nine array tasks were running: eight on `node-48-256g-8` and one on
`node-48-256g-14`. Twelve remained pending under the per-user job-count limit.
The eight node-8 processes each had a distinct physical core, IDs 16–23; their
allowed CPU lists contained the two siblings of that one core. Node-14's
process had core 0, allowed CPUs 0/24. No two inspected recovery tasks shared
a physical core. Later pending placements were not inspected in this snapshot.

Two raw counter snapshots about 10 seconds apart confirmed:

| Node | Processes | Sampled busy CPU clocks | CPU time over interval |
| --- | ---: | --- | --- |
| node-48-256g-8 | 8 | 2799.995 MHz at both reads | 10.65–10.66 s/process in 10.649160377681255 s |
| node-48-256g-14 | 1 | 3161.249 → 3179.557 MHz | 10.98 s in 10.977051686495543 s |

Small apparent utilization above 100% reflects sequential counter reads and
tick rounding. Existing thermal counters on these machines were nonzero, but
**none increased during the captured interval**. Package temperatures at the
second read were 61/40 °C on node 8 and 35/27 °C on node 14. The persistent
sub-GHz anomaly seen on the old node was absent in these recovery observations.

Actual second-attempt NanoShaper receipts also improved: qm36-La/Ca primary
meshes took 2.5341934263706207 / 2.537661850452423 s, compared with
36.952 / 36.233 s (rounded here) in the affected first attempts. Refined mesh
generation took 4.599–5.732 s. These are completed meshing measurements only;
they do not establish total solver speedup or numerical acceptance.

Raw capture files under
`workspaces/global_electrostatic_20260916/performance_observation/`:

- `node-48-256g-8_recovery_20260916.json`, SHA256
  `c24696a353751541ae8928886b5f2c8343f58c41a4cac2cde6000d82fb06f62e`.
- `node-48-256g-14_recovery_20260916.json`, SHA256
  `c98d4c1220ba8190394aecc9ca739358ca40681c4f6ad618f34b54ed3816116f`.

Per-task allocation records additionally reside under
`surface_recovery_array_v1/{00..08}/allocation_JOBID.json`. Continue the
unchanged recovery; no further technical intervention is indicated by these
observations.

## Grouped recovery placement

The twelve previously pending tasks now run as job1199974, using the existing
multiworker executor. Its admission receipt verifies twelve distinct physical
cores before starting. A subsequent read-only observation on node-48-256g-14
found thirteen active TABI processes on thirteen distinct physical cores
(the group plus one array task), at 2.800–3.094 GHz. No affinity, frequency,
priority or node-administration setting was changed.

Raw observation: `workspaces/global_electrostatic_20260916/resource_observation_group_v1.json`
SHA256: `28d504c011506aef5743c4a503be6a9044f3400124ff89aa908c4f33aef38c72`.
The twelve original pending array elements were cancelled before execution;
they incurred no solver attempt. The replacement waits for their termination
and retains the same scientific task manifests.
