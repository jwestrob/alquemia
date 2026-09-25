# Small LanM occupancy pilot — whole-node recovery, 2026-09-24

## Latest — 2026-09-25: two converged solvent cells, two vacuum recoveries queued

1216564 ended after5h52m21s: both ALPB cells passed collection; both vacuum
cells failed SCF convergence after500cycles. The independent reporter delivered
the actual2/4 result to the vault and local mail relay successfully. Read
`NATIVE_RESULT_1216564.md` and `DELIVERY_1216564.json`. These are numerical failures,
not scheduler or memory failures. All original artifacts remain intact.

Following Jacob's “two failed” / “proceed”, technical numerical recovery
**1217219** repeats only the two vacuum cells with their own converged same-metal
ALPB electronic states as initial guesses. Read [exact scope](SEEDED_VACUUM_PLAN.md).
The target Hamiltonian stays vacuum; atoms/coordinates/state/tolerances stay fixed.
Actual native restart markers must be observed. No additional solvent/MACE/DFT
calls or relaxed structures. Two workers use all allocated CPUs and node memory
with headroom. This is a separately named single-seed recovery, not proof of
a unique electronic state or a validated preference.

Runtime artifacts: `workspaces/lanm_global_occupancy_20260923/seeded_vacuum_v1/`.
Real-fixture preflight and missing-output collection passed in
`seeded_vacuum_preflight_v1` (no molecular execution there).
Agent `/root/lanm_completion_watch` watches start/native-restart/terminal status.
Independent reporterPID4160763 writes result/vault/email even without active chat;
receipt `DELIVERY_WATCHER_1217219.json`. No automatic further retry/expansion.

## Historical completed execution checkpoint

## Running and chat-independent delivery

**1216564 is actually RUNNING** on node-112-1500g-1: all112 CPUs are used by four
28-rank ORCA calculations, with1546754MiB node RAM and MaxCore10800MB/rank.
All four have entered actual native SCF iterations; no converged endpoints yet.
At10min26s, batch CPU time was17h49m17s and MaxRSS761154068KiB. This demonstrates
real CPU use; it does not establish convergence or scientific selectivity.

The preceding1216547 failed immediately because `SLURM_MEM_PER_NODE` is absent
for this scheduler's `--mem=0`. The fix reads actual `MinMemoryNode=0` and node
`RealMemory`; regression against real scheduler records passes. No chemistry
changed and that failed attempt ran no molecular evaluations.

User requested unattended execution through possible usage exhaustion. The batch
executes and collects independently of chat. A detached low-cost reporter is
running as PID3696728, recorded in `DELIVERY_WATCHER_1216564.json`; it polls only
this job, writes `NATIVE_RESULT_1216564.md`, a vault note and
`DELIVERY_1216564.json`, and emails the actual terminal result. Slurm END/FAIL
email is also enabled. Agent `/root/lanm_completion_watch` separately pings root
and writes `WATCH_1216564.json`; it does not duplicate email or submit work.
Agent messages can wait until the root session is active: do not claim that
these notifications can resume an assistant while its usage is exhausted.
Collection/reporting do not require an active assistant.

Only the four original native cells are running. No automatic resubmission,
new MACE/DFT or eight-source continuation. The invalid relaxation search remains
closed pending a physically valid motion policy. Inspect these actual receipts
on resumption before doing further work.

## Earlier recovery checkpoint

**Latest:** job1216461 failed after308seconds on64 allocated CPUs. The MPI slot
fix worked; ORCA then explicitly refused SCF because MaxCore2000MB was below its
5535.5MB vacuum /5735.8MB ALPB estimate. All four cells remain unavailable.
The PMIX warnings in output tails are not the identified fatal cause.

Jacob explicitly requested: “If we have the memory we should use more processes!
Just use what's available on the node!” He also requested completion watchers.
**Job1216547** now requests one exclusive memory-partition node and all its RAM.
At startup it divides ALL allocated CPUs among four MPI workers and derives
per-rank MaxCore from75% of allocated RAM (remaining25% covers other allocations).
The64-task scheduler request is an admission minimum; larger allocations use
their full CPU count:64→4×16,112→4×28,224→4×56. Each MPI launch fits the scheduler
slots and total simultaneous ranks equal allocated CPUs. No oversubscription or
scheduler/environment rewriting. All observed memory-node classes pass dry-run
resource checks. Input changes are MaxCore and MPI count only; source coordinates,
Hamiltonian, chemical states and SCF tolerances are identical. Cache keys change.

`native_feasibility_retry_v2/manifest.json` is created on the allocated node from
actual resources. `native_memory_preflight_v1` is dry-run only. The unchanged four
origin cells are the only molecular tasks; no new MACE/DFT or other sources.
Slurm END/FAIL emails are enabled; agent `/root/lanm_completion_watch` monitors
start and terminal status and pings root with actual outputs. Expected receipt:
`diagnostics/lanm_global_occupancy_20260923/WATCH_1216547.json`.

## Earlier first-run and first-retry checkpoint

Read [the first-run findings](FIRST_RUN.md). Job1213018 ran458seconds and failed;
its dependent1213040 was cancelled. Whole-protein MACE completed129 evaluations,
including a passing full/native force comparison, in451.4204568 worker seconds.
Both nonorigin proposals were rejected for stretched covalent C–C bonds.
There is no valid accommodated score or within-series preference result.

All four native origin tasks failed before energies because the original
one-task/32-thread allocation exposed too few MPI task slots. Technical recovery
**1216461** is submitted CPU-only with32 tasks,1CPU/task,200000MiB, four concurrent
eight-rank endpoints. `native_feasibility_retry_v1/manifest.json` has identical
input/coordinate hashes, scientific keys, states and implementation to the four
original tasks; only execution directories/concurrency changed. Dry-run and
syntax checks pass. No new MACE calls or DFT are included. Its collector retains
unavailable cells and explicitly leaves accommodation unavailable.

The other eight systems are not being run. Do not restart the old automatic
continuation: its completion gate does not require an admitted nonorigin geometry.
A physically valid motion policy must be established before expanding response
calculations. Do not loosen the covalent gate or count the large vacuum energy
drops as physical accommodation. PQQ production and the other session's PLM scan
are unaffected.

## Historical pre-execution checkpoint — superseded by the status above

Jacob explicitly requested a few two-ion and four-ion structures, beginning with
Hans-LanM, before extending across the library. He then requested email when a
meaningful result becomes available. The [frozen plan](PLAN.md) is the finite
execution scope. Existing PQQ production and manuscript delivery are unchanged.

## Prepared and checked

- Three complete chain-A structures: Hans8DQ2, Hans8FNR, Mex8FNS.
- EF1+EF2, EF2+EF3, EF1–EF4 for each source: nine conditional systems,18 La/Dy endpoints.
- No synthetic carve atoms, changed proton inventory, invented experimental labels,
  extra folding, new DFT or reserved SpyCI-LAMBS scoring.
- Five preparation tests and seven scorer tests pass on real pinned fixtures.
  These are preparation/parser/algebra checks, not successful molecular evaluations.
- One warm native MACE proposal per metal, all physical coordinates allowed within
  the declared small box. Every admitted geometry is cross-scored for both metals
  using native MACE + nativeGFN2(ALPB−vacuum).
- Dy physical multiplicities11/21 and native effective singlet are distinct,
  explicitly checked states. Multiple-spin coupling/ground-state validity remains
  untested. This is complete-monomer motion, not folding or dimerization equilibrium.
- Different source water inventories are retained. No raw cross-source energies
  or two/four-ion populations are compared.

## Jobs and exact continuation

At this checkpoint both jobs are pending; no new molecular evaluation has run.

1. **1213018**: first Hans8DQ2 EF12 proposal/common pool plus four native origin
   scalar calls. One H200,32 CPUs,200000MiB; at most24 native MPI ranks +8 MACE host
   threads simultaneously. Native/adapter full-source analytic forces are compared.
2. **1213040**, `afterok:1213018`, kill invalid dependency: verifies actual first
   complete MACE matrix, adapter equality, and four normal nativeSCF/state receipts.
   Only then runs the other eight systems on one H200/32 CPUs/200000MiB. It creates
   the exact remaining native scalar manifest and submits `run_native_matrix.sbatch`.
3. That future CPU-only stage uses32 CPUs/200000MiB, four concurrent eight-rank
   native endpoints, with four exact successful origin reuses. Maximum original
   matrix is54 composite endpoints/108 native scalar calls before deduplication;
   MACE optimization/qualification calls are counted separately. No retries or
   additional optimizations are automatically added.
4. CPU completion collects valid/failed cells, writes the conditional comparison
   and vault note, and emails the actual report to Jacob. Failed executor return
   codes do not suppress partial-result collection. Delivery is recorded as local
   relay acceptance, not assumed inbox delivery.

Do not duplicate these jobs or restart the old completed LanM pocket study.
Inspect actual receipts before interpreting results. A failed required cell/search
leaves that accommodated result unavailable; any valid static contrast remains a
separate field. No favorable-label gate controls the continuation or email.

## Artifacts

All paths below are beneath the repository root:

- `workspaces/lanm_global_occupancy_20260923/prepared_v1/manifest.json`
- `workspaces/lanm_global_occupancy_20260923/scoring_v2/manifest.json`
- `workspaces/lanm_global_occupancy_20260923/native_feasibility_v2/manifest.json`
- `workspaces/lanm_global_occupancy_20260923/AUTOMATION.json`

`scoring_v1`/`native_feasibility_v1` were preparation-only drafts, never submitted;
the v2 snapshot fixes missing-cell/search handling before any molecular call.
The native Hamiltonian and declared scientific search did not change.

Future result paths (not evidence of completed calculations yet):

- `workspaces/lanm_global_occupancy_20260923/RESULT_v2.json`
- `diagnostics/lanm_global_occupancy_20260923/PILOT_REPORT.md`
- `diagnostics/lanm_global_occupancy_20260923/pilot_result_email_receipt.json`
- Vault `agent-captures/2026-09-23_Nikasha-whole-chain-LanM-occupancy-pilot.md`

Root owns scoring/integration/jobs; Khoury's preparation and independent code
review are complete. Preparation committed73e9c7f. No remote push or baseline
promotion. Scheduler walltime defaults are not an agent-imposed compute budget.
