# Native whole-node execution recovery — 2026-09-24

Job1216461 entered ORCA STARTUP and GUESS, then all four tasks failed before SCF
with `Not enough memory available`. The configured MaxCore was2000MB; native
ORCA estimated5535.5MB for vacuum and5735.8MB for ALPB. Slurm records308seconds,
64 allocated CPUs (19712 allocated CPU-seconds) and batch MaxRSS222444248KiB.
This was an ORCA internal memory-setting failure, not evidence of SCF
nonconvergence or a scientific selectivity result. PMIX warnings also occur;
they do not supersede the explicit fatal error.

Jacob requested more processes wherever node memory permits, use of the whole
allocated node, and completion watchers. Recovery1216547 is CPU-only and requests
one exclusive memory-partition node with `--mem=0`. These are the documented
[Slurm whole-node CPU and memory options](https://slurm.schedmd.com/sbatch.html#OPT_exclusive).
The64-task admission minimum allows any current memory node; runtime reads the
actual granted CPU count, RAM and MPI slots. Four simultaneous workers divide all
CPUs equally. No environment variables describing the Slurm allocation are
rewritten, no oversubscription is enabled, and no other job is touched.

Actual node configurations checked before submission:

| Node class | CPUs | Allocatable MiB | MPI ranks per cell | MaxCore MB/rank |
|---|---:|---:|---:|---:|
| memory64, ordinary RAM |64|773914|16|9500|
| memory64, reduced registered RAM |64|677146|16|8300|
| memory112 |112|1546754|28|10800|
| memory224 |224|3094843|56|10800|

MaxCore uses75% of registered/allotted RAM across all ranks, rounded down to100MB.
The remaining25% covers allocations outside MaxCore. This is a resource setting,
not a project compute budget. Allocation/affinity/memory mismatches fail explicitly.
These values describe launch sizing, not measured scaling or achieved speedup.

The fresh manifest is prepared at job start under
`workspaces/lanm_global_occupancy_20260923/native_feasibility_retry_v2/`.
The four original coordinates, metal/charge/spin states, Hamiltonian and strict
SCF settings are preserved; MaxCore/MPI settings and cache keys change. The
fresh64-CPU dry-run under `native_memory_preflight_v1` passed, with exact
coordinate/input-body comparisons. Shell and Python syntax checks pass.
No additional molecular work was run in those checks.

The batch collects failed as well as successful execution and enables Slurm
END/FAIL email to Jacob. `/root/lanm_completion_watch` is separately watching
the job and reports start/terminal results to root, so a failed dependency
cannot silently suppress notification. No automatic additional molecular retry
or continuation is attached. Both earlier relaxed MACE geometries remain
rejected; four successful origin cells would establish native feasibility only.
