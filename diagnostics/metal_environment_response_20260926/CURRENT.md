# Field-aware metal response — restart checkpoint

## Active force-check batch — 2026-09-27

Root resumed authorized scientific work after actual completion wake was fixed.
Twenty fixed native analytic-gradient tasks are prepared/dry-run passed under
FORCE_CHECK_PLAN.md: two physical derivative modes, two signed step sizes,
repetition and joint rigid transform, both metals. Six real-artifact tests pass.
No new derivative result yet. Original saved gradients have net translation
residual below2.75e−6kcal/mol/Å; this sanity check does not replace the new tests.

**Owned job1219316** waits afterany PQQ1219220 (upstream1219217→18→19).
Exclusive344MPI task slots,20workers×17ranks=340, `--mem=0`; renderer uses full
scheduler node RAM with operational headroom. **Collector1219317** follows
afterany1219316. Workspace: workspaces/metal_environment_response_20260926/force_checks_v1/.
Manifest SHA2569c0a0c6f885c09f813c746c188e13cbed2373d6afdcf88e25991815bec745a9c.

**Wake watcherPID3158143 is armed**, event native-force-checks-v1, monitoring
both owned IDs and targeting root01a0a63a-ee36-7483-b726-1f7f5b7f75cc via the
tested `codex queue` interface. SUBMISSION.json retains exact launch/rearm
command and receipt path completion_event.json. Final output: FINAL_COLLECTION.json,
AUTO_REPORT.md and final_accounting.txt. On wake acknowledge the receipt, inspect
all force/refinement/rigid residuals, then decide whether expanded-region work
is justified. No automatic full hybrid/ML qualification or classifier promotion.

On restart inspect squeue/sacct and receipts; rearm watcher only using command
in SUBMISSION.json if needed. Do not relaunch successful chemistry. The previous
scout remains complete and immutable. No other-session job was modified.


## Latest terminal result — 2026-09-27

Scout1219207 and collector1219208 completed. All6native energies/analytic
gradient artifacts pass the existing parser; new directional derivative and
expanded-region tests remain pending. Thr159H rotation lowers electronic
energy for both metals: Ca -3.328669, La -2.975873 kcal/mol. Electronic response double
difference +0.352796 kcal/mol favors Ca along this specific perturbation; not an
affinity or classifier result. Native runtime276s,344allocatedCPUs=94944CPU-s,
0GPU-s; only4x16ranks were configured. Preserve actual cost and do not claim
whole-node utilization. Full hybrid and exact ML checkpoint remain blocked.

Read RESULT_1219207.json and reference_scout_v1/AUTO_REPORT_1219207.md.
No owned job remains active. Other session's new PQQ jobs1219217–1219220
are separate and must retain priority. Earlier queued entries below are dated
history. No new molecular submission was made during this status check.

Completion gap: the collector wrote files but did not notify/wake this chat.
Do not describe this as a completion-alert mechanism. For short pilots keep
an active monitoring turn through terminal collection; a detached collector
alone does not ensure the assistant resumes. No external message was sent.


2026-09-26. Root owns integration/execution. New user handoff authorizes the
staged Ca/La environmental-response scout and supersedes the previous discussion
pause. Historical whole-protein LanM vacuum subtraction remains closed.

## Current state (updated 2026-09-26 23:33 UTC)

- **Reference job1218751 queued afterany final discriminator collector1217591.**
  Six native analytic-gradient endpoints,64CPU task slots/256GiB, four16-rank
  workers, no GPU. All six cells currently unavailable/pending execution.
- **Collector1218752** runs afterany1218751, independently of this chat/login.
  It writes FINAL_COLLECTION.json, final_accounting.txt and AUTO_REPORT.md in
  the reference workspace. It launches no dependent scientific experiment.
- Upstream audit complete: code supports atomwise potential and field, but
  exact trained checkpoint and paper-compatible engine interface unavailable.
  Author question drafted in upstream/AUDIT.md and **not sent**.
- Real 1H4I preparation complete:54atoms,9087permanent charges, source-matched
  Ca−3/La−2 singlets; only Thr159HG1 rotates+10degrees in environmentB.
- Twelve unit/parser/allocation tests and five real preparation checks pass.
  No new molecular evaluations, native force qualification or predictive result.
- All three agents finished and released ownership; no duplicate job submission.
- Root owns interaction ledger, reference harness, numerical gates and submission.
- Production PQQ/PLM belongs to another session; preserve its jobs and dirty files.
- Reading anchor: git `87f1134`; existing unrelated tracked edits were present.

## Scheduling and interruption

Jacob says biotite shuts down Monday 28 September at 10 AM PST. Plan for the
earlier interpretation, 10 AM local Pacific, and do not rely on a surviving
login process. Research must queue **behind discriminator jobs**. Current
discriminator terminal collector is 1217591; other active chain IDs are
1217426,1217440,1217453,1217454,1217588,1217589,1217590. Requery queue and
dependency graph before submission; use explicit dependencies, not priority
changes to another session's jobs. Pending1218751 holds no running allocation.

Every submitted manifest, job ID, command, receipt and collection status will
be recorded here. No automatic rerun after restart: inspect scheduler and
receipts first. Keep partial and failed attempts. No email or remote push is
requested by this phase.

## Exact artifacts and next work

Workspace: `workspaces/metal_environment_response_20260926/reference_scout_v1/`.
Manifest SHA256: `c822b82da3857e45ead36a345d6130d8138833c56e065113dc83ed4ad3526696`.
Read SUBMISSION.json, COLLECTOR_SUBMISSION.json and DRY_RUN.json. Actual Python
dependencies and allocation-aware renderer are frozen under implementation/;
production scripts/environments were not altered. Native parser was exercised
on an actual archived embedded GGR output, not fabricated output.

Full additive QM/MM remains blocked by cross metal/PQQ repulsion/dispersion and
boundary-reference coverage. Scout measures only the declared **electronic
embedding component**, not total hybrid energy or affinity. PLAN.md freezes
its Hamiltonian, A/B rule and numerical tolerances. The complete model gate is
not passed by running a component test.

After1218751/1218752 finish, read actual collection/accounting. If complete,
declare finite native directional checks (MM plus mapped boundary), repeat and
rigid-transform cases, and expanded-region qualification before interpretation
or extension. No ML call without an exact accessible checkpoint.
No automatic4MAE/nonPQQ/LanM expansion, whole-protein optimization or promotion.
New molecular cost is zero while dependency-pending; preparation/code audit
time is separate and not a molecular timing benchmark.

Recovery commands (from repository root):

```bash
cat diagnostics/metal_environment_response_20260926/CURRENT.md
squeue -u jwestrob -o '%.18i %.50j %.10T %.30P %.30R'
git status --short --untracked-files=no
sacct -j 1218751,1218752 --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,MaxRSS,ExitCode
cat workspaces/metal_environment_response_20260926/reference_scout_v1/AUTO_REPORT.md
```

If the collector is interrupted, first inspect jobs/receipts, then collect to
a new filename without rerunning chemistry:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/metal_environment_response_20260926/reference_scout_v1/implementation/metal_environment_reference.py collect \
  --manifest workspaces/metal_environment_response_20260926/reference_scout_v1/manifest.json \
  --output workspaces/metal_environment_response_20260926/reference_scout_v1/RECOVERED_COLLECTION_1.json
```

If that filename exists, inspect it rather than overwrite it. Partial/failed
endpoints require an explicit fresh attempt, never an automatic rerun.
Startup exit75 means new discriminator work or Monday's shutdown cutoff
prevented molecular execution; read the job log. The cutoff is a shutdown
safeguard, not a project compute/time budget.


## Queue eligibility update — 2026-09-26

Jacob permits any available queue except test. Updated pending1218751 in place
to standard,memory,high-memory,standard-shared; afterany1217591 remains intact.
Suitable CPU nodes are idle; no GPU allocation is needed. No duplicate job or
scientific-input change. Actual scheduler update and before/after records:
`workspaces/metal_environment_response_20260926/reference_scout_v1/PARTITION_UPDATE_1.json`.
Original submitted batch script remains historical; this scheduler update is
the current partition policy. Collector1218752 remains standard-shared.


## 2026-09-27 recovery

Original1218751 exited75 before molecular execution because a new PQQ4A
batch was pending; collector1218752 recorded six missing cells. No chemistry
failed or ran, zero allocated core-seconds for the deferred scout. Queue now
empty; root resubmitted the unchanged frozen manifest as **1219207**, collector
**1219208**. No scientific input or historical collection overwritten.
See reference_scout_v1/RELAUNCH_20260927.json. New terminal files are
FINAL_COLLECTION_1219207.json, final_accounting_1219207.txt, AUTO_REPORT_1219207.md.
Original FINAL_COLLECTION.json/AUTO_REPORT.md remain the deferred attempt.


## Exclusive-node correction — 2026-09-27

Jacob reiterated: use full node RAM/threads, not the request as a memory cap.
Live1219207 acquired344CPUs exclusively but original request/manifest launched
4x16ranks with256GiB-derived MaxCore. All6tasks already started when inspected;
3Ca completed. Do not rewrite their frozen implementation or discard results.
Future renderer now verifies scheduler exclusivity and uses full RealMemory;
shared nodes remain allocation-bound. Future prepare requires explicit workers
and MPI ranks; execution checks all workers fit task slots. On this344CPU node,
6x57ranks uses342slots; memory is derived from8256990MiB with25%operational
headroom, not256GiB. Five resource tests pass, including actual captured scheduler
records in RESOURCE_POLICY_20260927.json. New run_reference_fullnode.sbatch
requests344slots/exclusive/--mem=0 for an appropriately prepared new manifest;
it is NOT a request to repeat the current six endpoints. Original submitted
script/manifest remain preserved. Policy is recorded prominently in AGENTS.md.


## Completion wake fixed — 2026-09-27

Scientific work was paused at Jacob's request while a subagent repaired actual
Slurm-to-root wake. Supported installed `codex queue` succeeded: direct probe,
then test jobs1219311(success)/1219312(intentionalfailure), automatically resumed
this root after final with no human nudge.90CPU-s total, no science. Read
diagnostics/slurm_wakeup_20260927/REPORT.md and COMMANDS.md. AGENTS now requires
an armed detached watch_and_queue.py monitor before yielding with owned jobs.
A saved report or subagent message is not enough. Watcher queues into the current
thread, then root consumes it after final; host restart requires rearming on
original job IDs, not redoing chemistry. Test gate passed, notification pause
is resolved. No new scientific calculation launched during this repair. Next
scientific work remains native MM/boundary directional-force qualification,
respecting current other-session PQQ priority and full exclusive-node resources.


## Scheduling correction — current ownership

Jacob clarified: do not raise priority over PQQ; **do not wait for PQQ**.
Cancelled only owned pending1219316/1219317 and replaced the submitted startup
guard with a normal-priority batch. Current worker **1219319**, collector
**1219320**, wake watcherPID3197471; event native-force-checks-normal-queue.
Same frozen20-cell manifest; no molecular work discarded. No other job or
priority modified. See force_checks_v1/NORMAL_QUEUE_SUBMISSION.json for exact
commands/rearm. completion_event_normal_queue.json is now the active receipt.
Earlier PQQ dependency instructions are superseded by this clarification.

## Native force result and grid diagnosis — 2026-09-27

Worker1219319 and collector1219320 completed; actual completion event automatically
woke root. All20 molecular cells complete;14/16 derivative/repeat/rigid checks pass.
Both rigid checks fail unchanged tolerances, while selected MM and boundary
finite differences pass. Read FORCE_CHECK_REPORT.md and RESULT_1219319.json.
A distinct eight-cell quadrature diagnosis is submitted as1219437, collector1219438,
8x43MPI=344slots, full exclusive RAM, normal priority with no PQQ dependency.
Frozen GRID_CHECK_PLAN.md; new grid_check_v1 manifest and four real-artifact tests
pass. Attempt to invoke pytest found it absent; unittest ran4tests successfully.
Wake watcherPID3344984 is armed, event native-grid-check-v1; exact rearm command
is grid_check_v1/SUBMISSION.json. On restart inspect receipts before submission.
Next: read grid_check_v1/FINAL_COLLECTION.json and actual numerical grid headers,
compare rigid residuals and response shifts. No automatic further refinement,
no full-hybrid/ML/biological validation claim. Original gate remains failed.

## Grid diagnosis complete — 2026-09-27

Worker1219437/collector1219438 terminal,8/8 endpoints complete. Wake event received
and ROOT_ACK.json written. Refined double difference0.352640263kcal/mol differs
by only-0.000155555; both rigid energy gates still fail. No further automatic grid
round, LanM optimization or expansion is launched. All owned jobs terminal.
Read GRID_CHECK_REPORT.md and updated REPORT.md; original qualification preserved.
Exact MACEPOL-EF weights/interface and full hybrid cross interactions remain blocked.
Next recovery command (read-only): python -c 'import json; print(json.load(open("workspaces/metal_environment_response_20260926/grid_check_v1/FINAL_COLLECTION.json"))["checks"])'.
Next scientific decision: specific residual diagnosis or obtain exact checkpoint;
do not resubmit completed34endpoints. Molecular workers total452016allocatedCPU-s,
zeroGPU; latest200552CPU-s. No external email or push.

## Continued native diagnosis — 2026-09-27

Jacob: "Proceed! Don't stop!" Continuing with a concrete SCF-stopping hypothesis,
not another grid sweep. Four refined A/rigid outputs retained density changes
far above nominal TightSCF density thresholds. New explicit all-criteria/forced
SCF protocol preserves grid, Hamiltonian, geometry, state and125iteration limit.
SCF_CHECK_PLAN.md frozen before execution;3new real-artifact tests pass,4old-grid
regression tests pass, dry-run passes. Manifest in scf_check_v1.
Worker1219450, collector1219451, watcherPID3695497; full344slots
(4x86MPI), full exclusive RAM, normal priority/no PQQ dependencies. Exact rearm
command in scf_check_v1/SUBMISSION.json. Read actual residuals and rigid checks
on completion; B response/strict finite differences remain unrun.

Expanded Thr159 region preparation is ready but unscored (63QM,9078MM), same
original boundary convention and physical state. Preparation agent owns its report
and vault note; no expanded molecular submission before interpreting native gate.
Public model branch/release recheck remains blocked on exact weights; no new
checkpoint or engine release found. See upstream/RECHECK_20260927.md.
Do not treat prior REPORT.md terminal checkpoint as current job state.
