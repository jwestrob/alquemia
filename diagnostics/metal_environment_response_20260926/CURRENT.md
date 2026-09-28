# Field-aware metal response — restart checkpoint

## Active diagnosis — 28 September 2026

Jacob approved explicit-f/frozen-f model diagnosis and ongoing scoped commits/pushes.
All campaign commits through776a13b are now on origin/main. No owned molecular jobs
are active; other-session PQQ jobs remain untouched. Root owns saved explicit-f
numerical diagnosis; finite_field_derivatives owns archived frozen-f evidence.
Read electronic_diagnosis_20260928/EXPLICIT_F_DIAGNOSIS.md and its pinned JSON.
Five real failed outputs audited. No unique wrong-occupation diagnosis; isolated
near-root stagnation differs from embedded gross divergence. Frozen-f force and
transfer qualification is not inherited from old converged singlepoints.
User now says shutdown tomorrow morning; exact scheduler shutdown time has not
been independently established. Maintain restart-ready records without assuming
additional compute availability. No automatic relaunch from older plans.


## PBE0 capability attempt stopped — 27 September 23:00 PDT

Root received the live near-root stagnation wake and stopped only worker1219868.
Twelve consecutive TRAH/NR macroiterations54–65 remained above5× the printed
1e-5 orbital-gradient tolerance, with an energy span4.83699e-7Eh. The residual
oscillates rather than improving reliably; nearly constant energy is insufficient
for a qualified analytic-force reference. No final energy/gradient is accepted.
The attempt ran about5hours before the alert; exact terminal cost follows in the
collector1219869 receipt. This was one isolated50atom explicit-f Dy sextet PBE0
calculation, not a protein relaxation or a preference comparison.

ROOT_STOP_DECISION.json in dy_pbe0_capability_v1 preserves the live evidence.
Keep collector and terminal watcher running. On completion record actual cost
and unavailable endpoint; do not restart this attempt or increase iterations.
Both native r2SCAN-3c and this PBE0 treatment encountered electronic convergence
problems in the isolated core. This rules out the protein field as a necessary
cause, but does not establish a unique cause, incorrect spin, or impossible Dy
chemistry. Numerical precision, electronic occupations and the source model remain
unresolved; the plateau alone does not distinguish them. No further method sweep
or dependent LanM expansion is justified by this result. PQQ remains untouched.

Next work uses saved outputs to define a specific electronic/numerical diagnosis
and a defensible reference strategy. The exact MACEPOL-EF weights remain unavailable;
there is no qualified field-aware ML candidate or within-series discriminator from
this campaign. Preserve the successful Ca/La environmental-response work separately.


## Current: isolated native Dy failed; one PBE0 capability test running — 27 September 18:02 PDT

The live-health event for 1219790 was reviewed and acted on. PModel explicitly
failed SCF convergence; HCore plateaued above the unchanged tolerance and root
cancelled the remaining work. Neither supplies an accepted energy or gradient.
The completed collector preserves PModel as invalid and HCore as missing its
final receipt; both were actually attempted. Worker 7010 s × 24 CPUs plus the
1 CPU-second collector cost 168241 allocated CPU-seconds, zero GPU time.
Read DY_SMALL_GUESS_RESULT.md. Do not restart this two-start experiment.

One separately declared PBE0-D4/def2-TZVPP explicit-f sextet capability endpoint
is now running: worker 1219868, collector 1219869, workspace
dy_pbe0_capability_v1. Startup verified: 24 MPI ranks on standard-shared,
--mem=0, observed available RAM policy; atomic PModel initialization is active.
This is not yet molecular convergence. Method and basis both differ from native
r2SCAN-3c. No preference result, matched La calculation or embedded continuation
is implied. Read DY_PBE0_CAPABILITY_PLAN.md.

Terminal watcher PID 1293984 and live-health v3 PID 1293986 are alive; exact
commands and pins are in SUBMISSION.json. Near-root stagnation now triggers a
review wake using the printed SCF tolerance. No automatic cancellation or
relaxed scientific acceptance. On wake inspect actual endpoint.out, electronic
state/gradient evidence and FINAL_COLLECTION.json; do not extend an unchanged
failure. All PQQ production, other agents' jobs and historical outputs remain intact.


## Active recovery: small Dy two-start diagnostic — 27 September15:55PDT

Jacob explicitly authorized proceeding before shutdown. Saved largeDy diagnosis
shows immediate molecularSCFdivergence, not failedatomicguess or observedbasis
linear-dependency. No unchangedlarge retry. New finite2cells use exactconsumed
50atomHansEF3Dy core, native r2SCAN3c/ECP28/sextet, isolatedHCore/PModel starts.
NoPAtom (archivedengineunsupported), newbasis/spin or biologicalscore.

Worker1219790/collector1219791:24sharedCPUs,2x12MPI,--mem=0,normalpriority.
workspace dy_small_guess_v3; exactSUBMISSION.json has BOTH terminalwatchPID3782325
and livehealthv2PID3782329 (earlyDIIS/TRAHpathology+partialcompletion alerts).
Noautomaticcancel or arbitrarytimecap; root responds to queuedhealthalerts.

Startup-only1219782/1219786 failed before molecularcalls; actualsharedSlurm omits
bothmemoryenvvars for--mem=0. Initialzero-envfix insufficient. Resourceprobe1219788
exposedexactcause; correctedprobe1219789 passed on actualnode. Newrenderer verifies
schedulerMinMemoryNode0, uses observedlocalMemAvailable bounded byRealMemory for
sharedmem0; exclusivefullnodepolicyunchanged. Allfailedattempts retained.
Sixresource tests and existing3state/parser tests pass, oneDyexecutiontest unrun.
Alloriginal195atomoutputs intact; other2sourcepreparations stay unscored.
ReadDY_SMALL_GUESS_PLAN.md and lanm_preparation_feasibility/DY_SCF_DIAGNOSIS.md.

## Terminal collection: Hans scout — 27 September

1219501 cancelled/1219502 completed; completionevent acknowledged. Actualcollector
retains2completeLa endpoints;2Dy attempted/unconverged lack finalreceipts and remain
unavailable. LaAsp85work−0.176654390kcal/mol; noLa/Dycontrast. Read
LANM_EF3_SCOUT_RESULT.md for evidence/cost/monitoringfailure. Total6,691,331allocated
CPU-s includingcollector. Noownedactivejobs, no retries/transfer submissions.
Next: saved-outputdiagnosis/modelselection; completion-onlywatcher is insufficient
for futurecostlySCF. PreservePQQ and allpreparedunscoredsources.

## STOPPED: Hans scout pathological Dy SCF — 27 September

At Jacob's runtime inquiry, job1219501 had run8h17m36s on224allocatedCPUs.
La_A/La_B terminated normally after2h47m/2h6m; both Dy SCFs remained pathological
with large TRAH residuals/negative orbital gaps. Root cancelled ONLY owned1219501
after preserving outputtails/reason in lanm_ef3_hans_scout_v1/EARLY_STOP_20260927.json.
Collector1219502 remains afterany; actual completionwake still armed. Preserve all
outputs and completed La receipts. No unchangedretry or transfer-source execution.
Other2source preparations are complete but unscored. Completion-only monitoring
missed prolonged nonconvergence; a runtime convergence/stagnation monitor is needed
before further costly reference work. No measured CPU utilization claim.

## Active: first explicit-state LanM EF3 scout — 27 September

Worker1219501/collector1219502, node-224-3t-2,224CPUs,4x56MPI,--mem=0. Actual completion
wakePID929432, event lanm-hans-ef3-explicit-state-scout-v1. Exact rearm command in
workspaces/metal_environment_response_20260926/lanm_ef3_hans_scout_v1/SUBMISSION.json.
Four endpoints only: La/Dy × source/+2degreeAsp85;195QM/1696field charges,
actualLa/La/Na spectators, total+6; new documented local peptide-boundary policy.
NativeLaECP46/DyECP28, physicalsinglet/sextet. No wholeprotein vacuum retry.
Seven preparation checks, independent review,7runner tests and actualdryrun pass.
State parser replays all6CaLa outputs and actualopen-shell format; Dyexecutionunrun.
Plan LANM_EF3_SCOUT_PLAN.md frozen before execution; explicitstate protocol
nikasha_lanm_ef3_embedded_native_explicit_state_v1. No inheritedbands/affinityclaim.
On wake inspect all4energies/gradients, actualspin/ECP/nativecomponents/SCF, projection
and costs before remaining2sources. Failures stay visible; no unchangedautomaticretry.
Savedcombined1H4I report complete: commonclassicalwork retainsdelta+.352795817;
electronicrigidfailure remains. Initial source/state audits are finished. Preparation agent now prepares the two
remaining declared sources (Hans8FNR/Mex8FNS), without molecular submissions;
execution expansion remains gated on this first scout.

Historical checkpoints follow. Latest entry above supersedes job status below.

## Read first: 27 September, classical matrix complete

All owned Slurm jobs are terminal. Three-core electronic response matrix completed;
classical component checks now pass200/200. No new classifier gain or full hybrid
qualification is claimed. Native rigid-rotation failure, absent solvent, missing
exact MACEPOL-EF weights, and newly diagnosed archived hydrogen/tail strain remain.
Read COMPONENT_CHECK_REPORT.md, THREE_CORE_REPORT.md and
hybrid_feasibility/component_checks/SOURCE_STRAIN_REPORT.md.

Active agent work: saved combined-force integration; consumed LanM EF3 preparation
audit; installed La/Dy electronic-state capability. No molecular jobs authorized
by these read-only tasks. Root owns future finite submissions and actual wake.
Do not rerun old scouts, serialize behind PQQ jobs, or revive whole-protein vacuum
GFN2. Entries below are chronological history; later completed records supersede
older queued descriptions. Latest completed job pair1219497/1219498.

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

## Expanded partition diagnostic queued — 2026-09-27

StrictSCF1219450/collector1219451 completed; automatic wake acknowledged.
Four outputs failed frozen summary-residual admission; all actually ended inTRAH.
Raw rigid error barely changes; close SCF remedy, no further grid/solver sweep.
See SCF_CHECK_REPORT.md and scf_check/FAILURE_DIAGNOSIS.md.

Independent coarse physical partition diagnostic now submitted as1219460,collector1219461,
watcherPID4023339, event native-expanded-partition-v1. Four real expanded Thr159
Ca/La A/B cells, original native DefGrid3/TightSCF;4x86MPI/full exclusive RAM,
normal priority/no PQQ dependencies. Tests5pass; dry-run passes. Immutable plan
PARTITION_DIAGNOSTIC_PLAN.md; manifest and exact rearm in partition_v1.
Known fine numerical failures remain, so this is not qualification, affinity,
optimization or classifier promotion. Small effects unresolved; source chemistry,
assembly/protons/waters preserved. New B has QM HG1 motion (not changedMMfield).
On completion read FINAL_COLLECTION.json: compare metal responses, partition
double-difference and HG1 projected forces, report actual costs and limitations.
No automatic new core/LanM relaxation follows. Production stays untouched.

## Expanded partition result — 2026-09-27

Worker1219460 and collector1219461 completed,4/4 actual endpoints; completion
wake acknowledged. PARTITION_REPORT.md/RESULT_1219460.json record delta small
+0.352795817 versus expanded+0.520767257kcal/mol, shift+0.167971439, same direction.
HG1 differential arc forces preserve direction at both endpoints; individual
components change with representation. Fine numerical qualification remains
unestablished, complete hybrid and exact ML model still unavailable. No affinity
or predictive gain claimed.452s/155488allocatedCPU-s; cumulative864128CPU-s,zeroGPU.

No owned Slurm jobs remain active. Root continues authorized work: preparation
agent audits matched4MAE/nonPQQ sources under original structural selection rule;
finite_field_derivatives reviews legacy cross-interaction parameters for a coherent
full hybrid path. These tasks are preparation/read-only, not molecular submissions.
Do not repeat completed numerical/grid/SCF matrices. Next decision comes from
transfer-preparation and hybrid-feasibility records, not an automatic optimizer.

## Second PQQ response diagnostic submitted — 2026-09-27

4MAE sixcell worker1219487,collector1219488,actual wakePID102822,event native-4MAE-transfer-v1.
Prepared exact matched80QM/8774PC canonicaldry source, Ca-3/La-2 physicalsinglets,
explicitLa368electrons with46ECP. Same structural rule selectsThr154HG1 PC2407;
only thisfieldcoordinate moves10degrees, fixedcore. Historical15P603 adduct/water
exclusions and longOH1.18042A remain explicit conditional-model limitations.
Four real-input checks and existingrunner dry-run pass.6x57MPI=342/344slots,
full exclusiveRAM/normalpriority/noPQQdependency. TRANSFER_4MAE_PLAN.md frozen;
transfer_4MAE_v1/manifest.json and SUBMISSION.json pin inputs/rearm.
Numericalqualification remainsunestablished; this is response transfer only,
not prediction, optimization or class improvement. Read terminal collection
and compare per-metal responses without absolute cross-protein totals.

Preparation agent continues1F6S exact atomiccharge/map inventory (no runs).
Hybrid feasibility found reusable metalLJ and exact-sourceproteinSystem; agent
exports unchanged proteinmechanics and audits exact-statePQQ GAFF2 crossLJonly.
No forcefield choice/completehybridenergy/optimization is silently promoted.

## Third response core queued — 2026-09-27

4MAE1219487/collector1219488 completed6/6; delta+.254914852kcal/mol, localresponse
not biologicalclass. Read TRANSFER_4MAE_REPORT.md; automaticwake acknowledged.
Exactsourceproteinmechanics andPQQcrossLJ nowprepared, completeadditiveledger
underreview; solvent-consistentrelaxation andexactML remainunsupported.

1F6S worker1219489/collector1219490,actualwakePID307513,event native-1F6S-transfer-v1.
Distinct normalized52atom core/1880ff19SBcharge state, Ca-1/La0, twoexplicitwaters
and11mappedcaps. Realpeptidebonds/NHverified. Sameproximalhydroxylrule selects
Thr86 at4.64009A O-core,+10degrees. Remoteperturb retained, no signalrescue.
Fourrealprepchecks andexisting6celldryrun pass.6x57MPI/fullnodeRAM,normalpriority.
See TRANSFER_1F6S_PLAN.md and transfer_1F6S_v1/SUBMISSION.json for exactrearm.
Oncompletion inspect FINAL_COLLECTION.json; complete3coredevelopmentreference
comparison with allnumerical/unavailableML/fullhybridlimitations visible. No
automaticlibraryexpansion orLanMrelaxation. Root owns sharedenergy/execution;
agentmechanics task onlypreparesinteractionledger, noenergycalls.

## Three-core reference complete; classical scout submitted — 2026-09-27

1F6S1219489/collector1219490 completed6/6; delta+.269010220kcal/mol. Actualwake
acknowledged. THREE_CORE_REPORT.md records all18primaryendpoints plus diagnostics;
cumulative54QMoutputs/1319928allocatedCPU-s,zeroGPU. These responses are notclass
decisions. Formalpointmetal retrospectivelycapturesroughlyhalf PQQresponse,closer
1F6S; see FORMAL_FIELD_COMPARATOR*. StaticpairedR regionpromotion shifts-1.071921A,
-1.239893B kcal/mol: responsecontinuity is not staticcontrastinvariance. This is
electroniccomponentonly, no calibratedscore or uniquephysicalcause claim.

Finite32configuration classicalcomponentmanifest prepared/dryrunpasses,4realartifact
noContexttests passed. SingleCa_A scoutworker1219495,collector1219496,actualwakePID473985.
Uses1sharedCPU because OpenMMReference intrinsicallysinglethread; no344CPUexclusive
reservation, --mem=0. Exactrearm in component_checks_v1/SCOUT_SUBMISSION.json.
COMPONENT_CHECK_PLAN.md freezes components/directions/tolerances before energies.
Read scoutreceipt/time/forces andpartialcollection; ifexecutable/admissible run
remaining31 tasks inparallelsharedallocation, reusingCa_A. No newQM/MD/optimization.
Fullsolvent/electronicnumerics/MLqualification remain unavailable. Ledger/classical
success cannot retroactively pass them. All agentimplementationtasks completed.

## Latest: classical matrix and source-strain diagnosis — 2026-09-27

Scout1219495 completed in11s/1CPU; finite components and305812KiB process peakRSS.
Large retained-bonded/LJ energies require source diagnosis, not relaxation. Agent
finds severe clash in experimentally missing/rebuilt A596 tail; detailed audit pending.
Frozen remaining31 classical configurations submitted: worker1219497,collector1219498,
32sharedCPUs,--mem=0,normalpriority; Ca_A cache retained. Actual wakePID568878,
event classical-component-matrix-v1, rearm command in component_checks_v1/MATRIX_SUBMISSION.json.
Inspect FINAL_COLLECTION_1219497.json and200checks on wake. No new quantum calls.
Root acknowledged scout event. Jacob requested continued work and campaign email;
email artifact/relay receipt saved in component workspace. Production untouched.

## Latest: classical derivative matrix complete — 2026-09-27

Worker1219497/collector1219498 terminal; actual wake acknowledged.32/32configs,
200/200frozenchecks pass. COMPONENT_CHECK_REPORT.md records exact residuals/cost.
Classical A/B response cancels between metals for this specific H motion; static
crossLJ Ca-minus-La is+0.382390829kcal/mol. Not a classifier or qualified fullhybrid.
Classical stage686allocatedCPU-s including scout/collectors; noGPU.
No owned active Slurm job. Preparation agent assembles saved combined-force evidence;
mechanics agent diagnoses missing-tail clash; wake-repair agent now audits installed
La/Dy reference state/basis capability. All tasks read-only/no molecular runs.
Native rigid failure and absent solvent remain explicit; no dry relaxation.

Final scheduler accounting: worker1219868 cancelled after18101s on24CPUs; collector1219869 completed in4s on1CPU. Total434428 allocatedCPU-seconds, zeroGPU. Collector retains endpoint unavailable after cancellation; no valid quantum energy or gradient.
