# One native GFN2 rank reproduces the energies at lower measured cost

Completed 2026-09-23. **One rank is the preferred execution candidate for these
tested compact contexts.** All 24 calls converged, with identical SCF iteration
counts and parameter exports across rank settings. One rank preserved every
component and paired contrast far inside the frozen tolerances while consuming
substantially fewer CPU resources. No production or active225 setting changed.

This is numerical/throughput qualification of **scalar** native GFN2 energies.
It does not qualify native gradients, alter chemistry, improve biological labels,
or repair the previously observed native force-consistency problem.

## Matched comparison

The same eight cells were executed at one, four and eight ranks: Q9Z4J7
(145 atoms) and Q89GY2 (184 atoms), actual precision-pilot adaptive_Ca geometries,
each metal Ca/La and each medium vacuum/ALPB. Exact source XYZ and electronic
input bytes were retained. Only runtime PAL changed; all jobs used the same
H200 host's CPUs, eight concurrent workers, no GPU, and the preserved one-thread,
no-binding runner. See the frozen [plan](PLAN.md).

| MPI ranks per cell | Calls complete | Largest cell difference vs fresh 8 ranks, kcal/mol | Largest paired R difference, kcal/mol | Gates |
|---:|---:|---:|---:|---|
| 1 | 8/8 | 1.24845e-9 | 6.40285e-10 | pass |
| 4 | 8/8 | 6.42057e-10 | 1.22237e-9 | pass |
| 8 | 8/8 | 0 (reference) | 0 (reference) | complete |

The predeclared gates were 0.1 kcal/mol **for every vacuum and ALPB cell** and
0.2 kcal/mol for each composite Ca-minus-La contrast. No cancellation hid a
failed component. Native MACE energies were exactly reused; changes in the
composite contrasts therefore come only from GFN2 transfer differences.
Fresh eight-rank values also reproduce every archived eight-rank energy exactly
at the parser's recorded precision. No energy was selected across rank settings.

SCF cycles are identical across all three settings: Q9 Ca vacuum/ALPB 29/19,
Q9 La 71/45; Q89 Ca 31/21, Q89 La 37/39. Actual exported native parameters,
charge/multiplicity/atom/electron accounting, normal termination, runtime PAL,
host and execution receipts were verified. No retries or failed molecular calls.

## Actual cost

| Ranks | Job | Allocated CPUs | Allocation wall, s | Allocated core-s | Reported CPU-s, distinct steps | Median cell wall, s |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1210743 | 8 | 18 | 144 | 42.520 | 10.296 |
| 4 | 1210744 | 32 | 43 | 1376 | 664.655 | 34.236 |
| 8 | 1210745 | 64 | 36 | 2304 | 1107.012 | 29.417 |

The one-rank batch took half the eight-rank allocation wall time and one-sixteenth
its allocated core-seconds **in this experiment**. Executor-only times were
13.069 / 39.596 / 31.381 s; summed actual task-rank-seconds were
80.061 / 1118.868 / 1811.754. These agree that fewer ranks were useful here.
They are separate from measured CPU usage and whole-allocation charges.

Each setting ran once, sequentially on a shared host. External load and startup
variation remain in the measurements; four ranks being slower than eight
illustrates why this is not a universal scaling curve. Allocation time includes
staging/collection; the final eight-rank job additionally generated the combined
comparison. Individual-cell and executor timings retain that distinction.

Total experiment cost: **3824 allocated core-seconds, 1814.187 reported CPU-seconds
across distinct job steps, zero requested GPU-seconds**. All job and step sacct
records were retained; parent and child accounting were not added twice. Host
memory requests were 16/64/128 GiB. Scheduler batch MaxRSS was 0/2455332/916168 K;
the first is unavailable as a meaningful peak measurement, and these step values
are not independent whole-allocation memory peaks. Local preparation, tests and
reporting CPU is additional and unmetered. Zero MACE/DFT/geometry calls.

## Implementation, verification and next step

A thin isolated adapter reuses the existing ORCA executor and parameter/state
audit. Each rank has a separate task/cache directory. Its wrapper collects
completed/failed statuses even if the executor exits nonzero; afterany dependencies
would allow the remaining predeclared settings to finish. No shared workflow,
default scorer, method reference or running225 manifest was edited.

Four [real-artifact tests pass](TESTS_final.txt), zero skips after execution:
exact rank-only recipe and source replay, paired coordinates/actual old receipts,
rejection of a corrupted real manifest, and actual final component/contrast
algebra. Initial preparation encountered an archive-schema difference between
own and cross-metal native cells; it was corrected to read the actual completed
pool matrix before any run directory or scientific call existed. Earlier failed
preparation/test logs remain preserved separately.

Recommend a deliberate **versioned one-rank execution option** for this compact
scalar stage, retaining eight-rank historical results and the existing physical
method. These two contexts support that numerical/resource choice; broader
sizes or gradients would need their own evidence. No rollout or further calls
were made by this branch.

- [Exact comparison](../../workspaces/native_gfn2_ranks_20260923/run_v1/COMPARISON.json)
- [Actual costs and accounting pins](../../workspaces/native_gfn2_ranks_20260923/run_v1/COSTS.json)
- [Finite design and task manifests](../../workspaces/native_gfn2_ranks_20260923/run_v1/design.json)
- [Compact artifact index](ARTIFACTS.json), [commands](COMMANDS.md).
