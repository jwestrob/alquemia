# Existing GPU evaluator agrees with the CPU on four real points

Launcher recovery5079 completed all four A0AC La-conditioned sample3 comparisons
on an NVIDIA RTX A5000. Every reported energy difference is exactly 0.0 kcal/mol;
maximum Cartesian force-component differences are 4.44e-14 to 7.11e-14 eV/Å.
All pass the frozen 0.01 kcal/mol and 0.001 eV/Å gates. These are two real
geometries, each evaluated as Ca and La, using the unchanged checkpoint and
native GPU evaluator. This qualifies this finite execution comparison, not a
new classifier or the as-yet unexecuted dynamic search adapter.

| Geometry / metal | GPU request seconds | Archived CPU request seconds |
|---|---:|---:|
| Fixed flip / Ca | 3.646060 | 15.300675 |
| Fixed flip / La | 1.171009 | 16.278494 |
| Relaxed La / Ca | 2.409549 | 16.477152 |
| Relaxed La / La | 1.443796 | 16.585711 |

The native GPU worker took 9.902658 seconds including 1.143243 seconds loading
the model, with four calls and zero failed requests. Peak CUDA allocation was
5,980,067,840 bytes and native host RSS 1,622,752 KiB. The complete scout,
including runtime/preflight overhead, took 46.235214 seconds internally; Slurm
records 46 seconds, one GPU and 12 allocated CPUs. Per-request timings above
are descriptive measurements, not a whole-workflow speedup. The native worker
uses one CPU thread; allocated CPUs are not measured CPU utilization.

Original5066 remains failed before any model call. Nested `srun` could not
confirm its batch allocation. V2 changes only process construction: it invokes
the pinned GPU server directly as a child of the granted one-GPU batch, with
matching parent/probe/worker cgroups and inherited device visibility. No SSH,
outside-allocation execution, CPU fallback, solver or model change occurs.
The underlying cluster reason for the nested-step rejection remains unresolved.

Authoritative artifacts are in
`workspaces/pqq_failure_causal_synthesis_20261015/asn_gpu_reuse_root_v2/`:
`CONFIG.json` SHA256
`e3a2a38b1ddbd18d15e471e8dd22e0e84fc1a252879bad92a3b6307f090a65d8`,
`execution/RESULT.json`, `execution/EXECUTION.json`, runtime/cgroup evidence,
raw native cell receipts and `gpu_summary_5079.json` (SHA256
`45c147d74cc9ccb58cc49376b8e00076d7a38fc20641a2f171f3a8301175a737`).
The agent's separate `asn_gpu_reuse_v2` binding is unsubmitted; do not run it.

Next: a minimal request/receipt bridge for the same evaluator and original
canonical search, retaining the completed CPU scout and A0AC as exact reuse.
The remaining23 CPU fallback stays unsubmitted to avoid duplicate chemistry.
