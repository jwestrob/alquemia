# CP2K memory recovery runner review

Reviewed 2026-10-23T21:20:01.348088+00:00.
**Source verdict: passes the bounded resource-only delta review.**
Configuration/wrapper/release binding remains a root-owned pre-launch check.
No runner, parser, native program or scheduler command was executed in this review.

## Exact authorities

- Recovery RUN.py: `9244ba6036bf220760e27bca95f1e4eb0e63e18e55e039600c84f3317961d7f8`.
- Original RUN.py: `03d9338b26e8bd507a9fd9bb4523ca6b8102e2597cebd7f582a67fc5f90bd0d7`.
- Original source review: `ca732271dc8819ad826e9368db3de0cd047e4232b3cce1fa4c1712f9e8d9d9fd`.
- Saved failed-attempt receipt: `83a2186691d477c79798ce5483ad97d6dea051090b41d46abfcff19df10ab288`.
- Original actual RUN_RESULT.json: `e2da1f0beef956cc53b8744982b5ee1db3d99f574c61a9db404d14ee7392ce1c`.

An exact text comparison, reversing only the three declared changes, restores
all original runner bytes: recovery docstring; native cap 7200 to 7073 seconds;
and actual `SLURM_MEM_PER_NODE == 262144` assertion. There is no other source delta.

## Findings

1. The memory request is correctly specified as per-node `--mem=256G`:
   262144 MiB. The runtime guard rejects absent/different per-node allocation
   before syntax or molecular invocation. It does not fabricate allocation or
   silently fall back to the old per-CPU default.
2. Exactly 24 allocated CPUs and 12 MPI ranks × 2 OpenMP threads remain required.
   Image, executable, command construction, environment, assets, input copying,
   parser, convergence criteria and molecular settings are unchanged by this
   source delta. The runner still contains one molecular parent invocation and
   no internal retry, method fallback, automatic extension or scheduler action.
3. Original job21618 failed OUT_OF_MEMORY with no admitted R/F/P/U. Its saved
   runner elapsed time is 126.61188986198977 seconds, while the native
   command itself used 118.65459037292749 seconds. The new cap satisfies
   `7200 - ceil(126.61188986198977) = 7073`; this conservatively subtracts the
   entire original controller elapsed time from the original molecular budget.
   Syntax checking, hashing and cleanup overhead are separately recorded; this
   is not a guarantee of a 7200-second total scheduler allocation including all
   support activity.
4. Inherited exact task-path identity, existing execution/receipt rejection,
   existing diagnostic-destination rejection, exclusive logs and atomic receipts
   preserve the old failed attempt. Timeout termination remains restricted to
   the newly created child process group. No cancellation or modification of
   unrelated jobs is present. Remaining quartet cells are not automatically
   launched by this runner and still require the pilot admission join.
5. The current AGENTS.md 23 October instruction explicitly requires
   standard_expand. Historical no-memory-override notes belong to the older
   OMol branch; they do not impose a new prohibition on this separately
   declared allocation correction.

## Pending concrete release joins

The fixed recovery configuration must pin this runner and the unchanged parser,
original decks, basis/ECP files, installed binary and image. Its task mapping and
result root must be fresh and disjoint from cp2k_quartet_v1. The shared wrapper
must request one node, 24 CPUs, `--qos=standard_expand`, `--mem=256G` and no
requeue, without extending scheduler walltime. Root must release only the
unfinished pilot at this stage and preserve job21618's failed evidence.
The source review supplies no new numerical/method or resource-sufficiency
claim; 256 GiB may still prove insufficient, and no further retry is automatic.
