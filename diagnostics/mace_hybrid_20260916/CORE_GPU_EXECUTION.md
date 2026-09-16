# Approved core checks on the standard GPU node

2026-09-16. Jacob: “let's try the core checks then, there's room on the gpu
partition.” This authorizes running the four already approved core tasks on the
older GPU node while retaining the full-protein work on H200.

Selected task IDs from the unchanged scientific manifest:

- `1h4i_qm33_La` and `1h4i_qm33_Ca`: 47 atoms each.
- `1h4i_qm36_La` and `1h4i_qm36_Ca`: 54 atoms each.

Request one RTX A5000 on `node-128-512g-8gpu-1`, partition `gpu`, 16 CPUs and
64474 MiB host RAM (one eighth of configured node memory, rounded to MiB).
The same installed checkpoint, float64, realspace electrostatics, coordinates,
charge/spin and energy/force operations apply. This is scheduling the existing
four tasks, not adding scientific evaluations or a new model variant.

The H200 job must wait for the core job to end, because the existing executor
locks the shared campaign. It then reuses verified core results and performs
the remaining eight full-protein calls. Failed attempts remain visible;
hardware-specific execution receipts distinguish the two allocations.
Baseline/defaults remain unchanged. No new DFT or new biological cases.

## Completed execution

1200302 failed during unused reciprocal-grid construction;1200306 failed on
missing checkpoint metadata. Both failures remain archived. Versioned repair
v3 completed all four cores in1200308 on the prescribed A5000 allocation.
[Results](CORE_RESULTS.md) include all costs and the failed partition diagnostic.
1200207 was replaced by1200309 with the repaired implementation and completed
core cache; same approved full-system tasks, H200 allocation and memory policy.
