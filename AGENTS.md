# Agent entry point

## Current Pro review — 2026-09-26

Read [the current model-selection brief](docs/PRO_REVIEW_LANM_20260926.md).
Both seeded vacuum retries1217219 failed after500iterations; all owned LanM
calculations are terminal. The whole-protein MACE search also yielded no admitted
relaxed geometry. Jacob wants scientific model selection with Pro before more
engineering. PQQ execution belongs to another session. No new LanM run is queued
here. The restart/queued checkpoints below are historical; do not relaunch them.

## Restart recovery — 2026-09-25

Start with [the compact restart handoff](diagnostics/lanm_global_occupancy_20260923/RESTART_HANDOFF.md).
It records the actual results, queued1217219, approved two-cell numerical recovery,
exact resume/watcher commands, concurrent PLM ownership and vault location.
A host reboot can kill detached reporters and agent sessions even if Slurm jobs
survive. Inspect the scheduler and receipts before resubmitting anything.
Vault copy: `agent-captures/2026-09-25_Nikasha-restart-handoff.md`.

## Active small whole-protein LanM occupancy pilot — 2026-09-23

Jacob authorized a few two-ion/four-ion structures, starting with Hans-LanM,
before broader within-series work. Read
[current state](diagnostics/lanm_global_occupancy_20260923/CURRENT.md) and
[frozen scope](diagnostics/lanm_global_occupancy_20260923/PLAN.md).
Job1213018 executed on September24: whole-protein MACE passed native force
qualification, but both nonorigin proposals failed covalent geometry checks.
All four native origin tasks failed at MPI startup; dependent1213040 was cancelled.
CPU recovery1216461 also failed: ORCA's MaxCore2000MB was below its measured
SCF requirement. Jacob requested full use of allocated CPUs/RAM and watchers.
Recovery1216547 failed before any molecular work (absent optional memory envvar).
Fixed1216564 completed with two converged ALPB cells and two vacuum SCF failures
after500cycles. Technical recovery1217219 uses each metal's successful ALPB state
to initialize its unchanged vacuum target, two tasks sharing the full allocation.
Read SEEDED_VACUUM_PLAN.md in the diagnostic directory. No solvent cells repeat.
Agent `/root/lanm_completion_watch` pings root; chat-independent reporterPID4160763
collects, writes result/vault and emails terminal results. See current state for
receipt paths. Slurm END/FAIL is also enabled. User may exhaust assistant usage;
compute and reporting continue independently. No new MACE/DFT or eight-system
continuation runs. Do not assume a pending agent notification wakes an idle root.
Read [first-run findings](diagnostics/lanm_global_occupancy_20260923/FIRST_RUN.md).
Do not restart the old continuation: its gate does not require an admitted
nonorigin geometry. Accommodation and within-series preference remain unavailable.
PQQ production/SOP and all prior outputs stay intact.
Root owns scoring/execution; Khoury's preparation/review is complete. No reserved
SpyCI-LAMBS outcomes are opened and no library campaign is launched.

## PLM operations packaged — 2026-09-23

Jacob requested the final PQQ SOP and manuscript draft, with Spicy-Lams inventory
in parallel. Read [the SOP](docs/PLM_PQQ_SOP.md) and
[packaging report](diagnostics/plm_pqq_delivery_20260923/REPORT.md).
Fresh preparation of all six consumed PLM sources exactly reproduces the archived
envelope state/coordinates/maps. The operations package reuses existing science;
no new endpoint scoring, full-cohort submission or default promotion occurred.
Jacob plans to start the scan. Preserve the current production/static and explicit
DFT access, individual failures, and the separate recovery evidence.

The [Spicy-Lams inventory](diagnostics/spicy_lams_inventory_20260923/REPORT.md)
and vault note locate 616 assay sequences and full mature-chain models for all
16 reserved-panel proteins. No new within-series scoring or reserved-outcome
unblinding occurred. Whole-protein response is
the intended next research direction, but occupancy/assembly must be established
per construct and conditions rather than imposing four ions on every LanM.
The older overnight jobs below are complete; do not restart them.

## Completed overnight discrimination work — 2026-09-23

Jacob authorized discretionary contained experiments and parallel agents. The
round is complete; all owned jobs and collectors are terminal. Read
[the delivery](diagnostics/overnight_discrimination_20260923/DELIVERY.md) and
[checkpoint](diagnostics/overnight_discrimination_20260923/CURRENT.md) before
new work. Preserve production/static, explicit DFT and all original artifacts.

Strongest completed ten-fold-context/adaptive candidate: 207 correct, zero wrong,
1 inconclusive, 17 unavailable among 225 consumed structures; both released
errors repaired. Same-solver static pockets give 199/1/8/17, demonstrating a
contribution from accommodation. All 94 supported three-source subsets remain
correct. [Full result](diagnostics/strict_native_transfer_20260923/REPORT.md).

The practical 4.3 Å three-source envelope has now completed its full100 test:
primary 91 correct/9 unavailable, including one actual failed scalar cell that
blocks three triples. Two prescribed same-metal electronic restarts converge
and agree; the separately named recovery sensitivity restores 94 correct/6 old
exclusions. The original primary result remains unchanged. Four of104 individual
source/context pairs are inconclusive, versus one under strict tenfold; the
separate A0A3Ca3 pilot also becomes inconclusive. Fold spread does not uniformly
shrink. Keep these limitations alongside the complete median-call fidelity.
[Primary](diagnostics/motion_envelope_transfer_20260923/REPORT.md),
[recovery](diagnostics/native_failed_cell_recovery_20260923/POOLED_SENSITIVITY.md).

The actual two-PLM-triple integration completes all6 sources in251 allocation
seconds on1 H200/32 CPU, with183 MACE/72 strictGFN calls. Both medians remain
Ca-supported on the candidate scale; their true preferences are unknown.
Short-contact relief reduces the two source-score ranges substantially.
[Physical response](diagnostics/plm_envelope_response_20260923/REPORT.md).
The [reusable opt-in API](diagnostics/pqq_three_source_envelope_api_20260923/COMMANDS.md)
accepts explicit complete prepared triples. Interfacev2 preflights are unexecuted;
actual integration remains v1. Automatic restart is not part of the API.
Root's [overlay](diagnostics/plm_candidate_overlay_20260923/REPORT.md) appends only
these two actual groups to176 unchanged original protein/DFT/expression records.

Six-angle expansion failed to add useful discrimination and is closed. Strict
solver repetition on all old DFT-tested donor motions changes response by at
most0.081 kcal/mol; solvent-added DFT disagreement remains. LanM series transfer,
standalone-solvent replacement, generic scaffold proposals and solvent-gradient
continuation are also closed. Do not rerun them from older checkpoints.

All three agents have completed reports/tests/vault notes and released ownership.
Root owns delivery, shared docs, the joined tables and substantive email (sent
with the physical-response figure). No production promotion, cohort rescore,
new DFT/folding, remote push or goal-status change occurred. Older text below
is historical and does not authorize restarting completed campaigns.

## Restart/scaffold round completed; LanM follow-up requested — 2026-09-22

Read [the completed round](diagnostics/scaffold_restart_round_20260922/REPORT.md)
and [current ownership](diagnostics/scaffold_restart_round_20260922/CURRENT.md).
All seven jobs/collectors are terminal. Native continuation removes the Q88La
4.805 kcal artifact; no full predictive restart comparison has run. All80 saved
seed pairs for the proposed four-pool follow-up exist; that follow-up is not
automatically authorized by the completed diagnostic.

Corrected collective scaffold proposals score8/8, but transferred adaptive calls
degrade7correct/1wrong→3/5 and both A0AC La folds fall below two Ca controls in raw
ordering. Reduced fold spread does not rescue this result. Close this generic
protein-only proposal rule; do not enlarge its bounds or refit on these eight.
The reusable Cartesian shared-pool primitive remains available. Production stays
unchanged. Preserve failed CUDA, false-stereocenter and seed-only attempts.

Jacob then requested a subagent for within-lanthanide LanM discrimination.
Water_basins owns `diagnostics/lanm_series_followup_20260923/`: recover actual
last tests, source/label compatibility and electronic-state support. Root asked
whether to prioritize La/Dy solution affinity or La/Lu resin selectivity; these
are distinct evidence types. No new series calculations have launched. Do not
inherit La/Ca bands, silently force open-shell Dy into singlet adapters, or call
raw cross-element core energies affinity. Root owns integration; other branches
are complete. The broader improvement goal remains open.

## Completed parallel pilots — 2026-09-22

Jacob approved alternative starts, solvent-guided proposals and direct local
basin-width tests. All executions/collectors are terminal; read
[the final report](diagnostics/nikasha_parallel_pilots_20260922/REPORT.md) and
[checkpoint](diagnostics/nikasha_parallel_pilots_20260922/CURRENT.md). Neither
search policy adds discrimination. Two of four curves pass integration gates,
with small width contributions, but electronic-solver continuity is unqualified.
A reproducible4.805kcal GFN2 jump occurs at almost identical Q88 geometries;
large solvent-component jumps also appear on finite basin grids. No remedy or
new classifier gain is claimed. Production stays unchanged. Restart, standalone
xTB, collective-scaffold and redox follow-ons remain proposed-only. Do not rerun
these completed manifests or launch those separate proposals from this checkpoint.

## Latest completed phase — 2026-09-22

Jacob: “lets address all that, including adaptive accommodation.” Read
[the three-branch plan](diagnostics/nikasha_next_phase_20260922/PLAN.md).
All three branches and their jobs/collectors are finished. Read the
[final result](diagnostics/nikasha_next_phase_20260922/REPORT.md) and compact
`RESULT.json`. Adaptive accommodation now has a complete canonical reference and
full225 transfer: on204 matched sources, released200correct/2wrong/2inconclusive
becomes202/1/1. It fixes the A0ACD6B9F2 error and two inconclusives, retains all111
available Ca-class calls, but loses three sources versus released coverage and
adds cost. Pursue this branch; no production promotion. Consistent context alone
fixes a different error but adds six abstentions; CPCM yields no converged metal
pair. The proposed combined-context/adaptive test is not executed or authorized
by this completed phase; agree its scope before launching.
Do not rerun these manifests or collectors. Production fastMACE/GFN2 and explicit
DFT remain intact. No new DFT, folds or PLM rescore occurred. Root integrated the
three agents' work; their tasks are complete. Older checkpoints are historical.

## Latest completed delivery — 2026-09-22

Current name **Nikasha**; released fastPQQ, explicit DFT and historical commands
remain intact. Read [delivery](diagnostics/nikasha_recovery_20260922/DELIVERY.md).
All September20 preservedDFT/proposals and September22 common-pool/adaptive/joint
jobs are terminal. No active collectors or new submissions are needed. Independent
CC failed without a complete metal pair; do not retry it. The separately gated
candidateDFT campaign remains **dry-run only**.

The released static context composite remains strongest on the completed matched
structural challenge:203/207 correct versus DFT198/207, from25 consumed groups.
Both retain all94 available La-conditioned triples; these are correlated subsets.
Shared pooling did not earn routine use. Full225 report and frozen reference:
[shared-pool result](diagnostics/nikasha_shared_pool_20260922/FULL225_RESULT.md).

Four-angle adaptive scoring gives29/30 available. On24 matched canonical members,
static24correct versus adaptive22correct/1wrong/1inconclusive; missingMMOL1770
prevents its calibration. Jointmetal/four-angle pilot1209901–03 gives7/8 valid
candidates,3/4 pools, real native-force relief and no changed available class calls.
Both missing candidates are oversized optimizer-trial failures, not bad protein
chemistry. See [angular](diagnostics/adaptive_accommodation_20260922/CANONICAL_POOL_COMPARISON.md)
and [joint](diagnostics/adaptive_metal_20260922/FINAL_RESULT.md) reports. No promotion.

Scaffold feasibility is complete: exact standard-protein parent maps/parameters
exist, but a defensible capped-local subtraction and exterior metal/PQQ coupling
do not. No scaffold energies/corrections were manufactured or run. Its optional
parent diagnostic is proposed only; do not treat its presence as a queued job.

PLM export retains176proteins/200genes, originalDFT identities and actual transcript
measures. Two editable scientific figures and methods outline are ready under
diagnostics/nikasha_manuscript_20260922. No MACE cohort rescore, remote push or
publication occurred. All agents' bounded tasks are complete. The larger goal
of better discrimination remains open; more motion alone is not demonstrated
utility. Recover actual results before proposing another model/version.

## Latest: accommodation mechanism supported; test classifier robustness — 2026-09-20

Original native donor validation is complete: both PLM differential responses are
confirmed. Full composite optimization1204162 is also complete, but zero pairs
pass its frozen minimum/curvature gates; its native follow-on is **dry-run only**.
The cheaper MACE-proposal/actual-composite-selection experiment1204169/1204171
completed all30sources. Reference classes remain ordered; one Ca reference becomes
inconclusive on old bands. This is engineering progress, not yet an accuracy gain.
Read [proposal result](diagnostics/accommodation_nonlinear_20260920/PROPOSAL_REPORT.md).

The next [declared challenge](diagnostics/accommodation_nonlinear_20260920/FOLD_TRANSFER_PLAN.md)
uses every225primary noncanonical reference fold, with the same physical rule and
all unavailable inputs retained. A distinct canonical-only development reference
was frozen inbdec940; old-band transfer remains explicit. All208 supported physical maps are ready. Reviewed GPU1204185 now runs
the unchanged proposal rule (416tasks,415available origins). Khoury may launch its
four finite solvent shards after actual dry-runs; no further permission is needed.
Adapter7809fe0/review32ef387 and seven real-artifact tests passed. Water_basins is
finished; Khoury/root own execution and results.
Second_shell owns final existingDFT comparisons after1203976–1203979 complete.
All production defaults and PLM results remain unchanged; preserve independentCC.
Read [current checkpoint](diagnostics/accommodation_goal_20260920/CURRENT.md) before
submitting work, and do not rerun completed experiments.

## Latest: three-fold report works; physical response still under test — 2026-09-20

Optional `affordable_workflow.py standard ensemble` is committed as11d9669;
read [its commands and limits](diagnostics/pqq_ensemble_20260920/REPORT.md).
Every available consumed-reference triple is correct with context+solvent, and
the simpler native core matches that result. The two isolated PLM development
cases remain Ca-like across all three folds; [report](diagnostics/plm_fold_sensitivity_20260920/REPORT.md),
commit82ed189. Their biological labels remain unknown. No ensemble default or
production PLM rescore occurred.

NativeDFT1203771 confirms the small control response and the first large PLM
calcium response; lanthanum/differential validation remains incomplete. The four
preserved-baseline jobs1203976–1203979 remain active. Both have automatic collectors;
do not duplicate them. The same100-triple DFT comparison is implemented in11a41b0;
run its [final command](diagnostics/accommodation_fold_DFT_20260920/THREE_FOLD_DFT_REPORT.md)
after final_collection.json exists. Root owns the current checkpoint; all three
agents' original scoped tasks are completed or finishing reports. Reactivate
water_basins for full native response and khoury_benchmark for final matched DFT
comparison when actual collections exist. Preserve independentCC1202429.

## Current: fast PQQ released; structural robustness under test — 2026-09-20

The authorized fast PQQ source-to-score release passed and is committed as
8dd81ba. `affordable_workflow.py standard` chooses the named fast mode for
compatible explicit PQQ source requests; explicit `dft-reference` and existing
water-prepared baseline remain available. Read
[release report/limits](diagnostics/pqq_fast_release_20260920/REPORT.md) and
[commands](diagnostics/pqq_fast_release_20260920/COMMANDS.md). This release proves
source reconstruction and retained reference fidelity, not broad accuracy gains.

Root is testing all250 saved, known-reference AF3 samples across both folding
metals under [frozen plan](diagnostics/accommodation_goal_20260920/FOLD_ROBUSTNESS_PLAN.md).
Primary comparisons exclude the25 original calibration geometries. Root owns
`accommodation_folds.py` and execution; second_shell owns the comparison/report;
khoury_benchmark audits source/H-preparation failures. Water_basins owns the
separate solvent-consistent physical donor profiles and preselected DFT checks.
Do not change bands or infer PLM labels from predicted XoxF names. See the
[current checkpoint](diagnostics/accommodation_goal_20260920/CURRENT.md) before
submitting work; preserve CC1202429 and all independent jobs. Jacob requests
updates led by motivation and practical benefits, with detailed numbers in reports.

## Active overnight goal — 2026-09-20

Jacob authorized autonomous pursuit of improved metal-accommodation discrimination:
"Set a goal ... run on this overnight ... If you see promising new directions,
pursue them." Read [goal/ownership](diagnostics/accommodation_goal_20260920/GOAL.md).
The fast-PQQ source-to-score release is separately authorized for conditional
compatible-PQQ promotion; second_shell owns that integration. Root studies physical
donor organization; water_basins owns the solvent-consistent gradient pilot;
khoury_benchmark owns matched source controls. PLM results are predictions, not
labels. Existing PLM geometry screening already identified compressed sites; reuse
that work rather than rediscovering it. Speed or preparation alone is not goal
completion. Preserve all unrelated edits, CC1202429 and its completion monitor.

## Latest: usable opt-in composite scorer — 2026-09-20

Jacob said "pursue!". Numerical tightening passed all6contexts:24nativeGFN2
endpoints, maximum score shift0.025545kcal below unchanged0.20tolerance. Ordinary
SCF remains unreliable even after actual orbital restart; no complete cross-solver
pair and no broad qualification claim. Read [numerical result](diagnostics/compact_qualification_20260920/REPORT.md).

The fresh four-site prepared-input MACE+solvent path completed in117s on1H200/
32CPUs/200000MiB, about27–31s per site; all scores repeat within6.4e-10model-kcal.
Use [scanner commands](diagnostics/compact_scanner_20260920/COMMANDS.md) and
`scripts/compact_solvation_scanner.py` for opt-in supported prepared pairs. Old
MACE energies are optional. PQQ bands require compatible PQQ context provenance;
generic sites retain raw contrasts. Current commands use the corrected32-task
MPI layout; the first layout attempt failed before GFN2 SCF and is preserved.
Total continuation cost including all failures:21,248allocatedCPU-s/160GPU-s.
Default DFT/water-preparation workflow and calibration remain unchanged. Do not
rerun completed qualification/timing pilots. CC1202429/monitor remain independent.

## Latest research result: compact MACE + solvent transfer — 2026-09-20

Completed approved native GFN2 ALPB-minus-vacuum correction on the fixed33-case
core/context inventory:264/264 primary endpoints converged. Core native MACE
alpha/GGR improves2/6→6/6; both representations retain25/25PQQ plus3/3consumed
crystal calls. Context PQQ gap improves2.315→5.076model-kcal. Standalone GFN2 does
not achieve the same ordering. This is useful development evidence, not broader
accuracy than the promoted DFT baseline. Eight alternate-SCF checks failed;
numerical qualification remains unavailable. No new DFT/MACE calls or default
change. Read [result](diagnostics/compact_solvation_20260920/REPORT.md) and
[commands](diagnostics/compact_solvation_20260920/COMMANDS.md); do not rerun the
completed panel. All attempts cost37,376allocatedCPU-s, zero newGPU-s. Pursue the
challenger, resolve numerical reproducibility before promotion. Charged-PQQ DFT
1202478 is also complete; CC1202429 remains separate under its existing monitor.

## Current: water preparation promoted — 2026-09-19

Jacob authorized promotion and discretionary next experiments: “go ahead and
promote it ... make sure you have an eye towards improving the classifier.”
Use `scripts/affordable_workflow.py baseline` for the new versioned prepared-site
workflow, default `contextual_if_supported`; explicit `original` remains available.
Read [promotion report](diagnostics/water_promotion_20260919/REPORT.md) and
[executable commands](diagnostics/water_promotion_20260919/COMMANDS.md).
Supported wet amide-v3 sites receive MACE contextual water-H preparation before
unchanged native DFT. Original/prepared scores are retained side by side. Dry PQQ
is exact identity with its own released bands. Unsupported wet cofactors/legacy
cores fail explicitly; raw legacy inbox watchers have not been migrated/restarted.
No occupancy, entropy, environmental scorer or new absolute threshold is promoted.

Parallel accuracy findings and jobs are in
[current status](diagnostics/accuracy_tracks_20260919/STATUS.md).
Static complete second-shell context improves every alpha/GGR structural margin;
compact native OMOL retains25/25PQQ plus3/3consumed crystal transfers, with a much
narrower calibration gap. No context-scorer promotion. Source-supported geometry
proposals did not improve these margins; sampled internal proton transfer is
uphill. Independent parvalbumin water preparation completed, but its La/Ca
labels remain unresolved. Coupled-cluster and charged-PQQ native checks are
separate ongoing diagnostics; inspect live receipts before doing anything.

The user explicitly authorized parallel chemical-state, environment and
independent electronic investigations. Metal specificity includes noncatalytic
binders such as LanM. Catalytic competence and imitating unchanged DFT with a
faster model are not the accuracy target. Preserve active jobs and other agents'
files. Older entries below are dated development checkpoints, not instructions
to rerun completed work or undo the authorized promotion.

## Active: improve the MACE-assisted scanner — 2026-09-19

Jacob's active goal is in [parallel goal](diagnostics/parallel_discriminator_20260919/GOAL.md).
He owns remaining PLM folds; prioritize classifier utility before the provisional
Tuesday shutdown. Context-prepared masked MACE improves alpha/GGR 2/6→6/6
structural comparisons, with unchanged dry PQQ inputs; one biological comparison.
Read [scanner result](diagnostics/hydration_scanner_20260919/REPORT.md).
MACE donor response plus DFT also preserves 25/25 PQQ and improves consumed direct
ordering without fitting to those labels; see diagnostics/response_probe_20260919.
Second-shell and wider water-basin jobs and continuations are active; inspect the
live queue and their diagnostics. Preserve existing jobs and other agents' files.
These are progress toward the goal, not completion or default promotion.

## Current: coupled water response complete — 2026-09-19

27native DFT checks and all MACE tasks completed. All11trial steps lower native
energy and pass energy checks;3/4endpoints meet water-coordinate minimum criteria.
1F6SCa remains nonstationary (projected max0.70245 versus0.4 limit) and retains
its failed soft-mode check.7/8initial direction checks pass; all4response paths pass.
This supports cheap water-relaxation proposals, not harmonic entropy or occupancy.
Thermal model extents exceed the validated local domain; missing free-energy terms
remain null. Baseline/default/PQQ unchanged.71real-fixture tests pass, zero skips.

Read [coupled report](diagnostics/hydration_basin_20260919/REPORT.md),
[actual tables](diagnostics/hydration_basin_20260919/export_final_v1/TABLES.md), and
[commands](diagnostics/hydration_basin_20260919/COMMANDS.md). The finite continuation is
complete; do not launch a third recentering round or rerun completed jobs. One
qualified6IP9La endpoint was explicitly reused. Exact new allocation cost:
178496core-s/840GPU-s. Final report uses frozen analysis_implementation_v4.
Four initial collection equality failures were recovered without new chemistry.
The next research issue is wider water basins and complete state accounting,
not another identical local step or an invented entropy offset.

Older native orientation comparator1201825 for6IP9 remains live; preserve it and
collect it when finished.1201824 is complete with recorded geometry-tolerance
failures. These older jobs are separate from the completed coupled-response pilot.
Prior radial/occupancy/preparation results remain preserved; older summaries below
are historical checkpoints. Current work and vault note are recorded in SESSIONS.

Read [docs/AGENT_PIPELINE.md](docs/AGENT_PIPELINE.md) for current operations,
protocol choices, interpretation limits and runnable examples. Older root
documents retain dated material; the current guide identifies what supersedes it.

Before acting, inspect git status, [SESSIONS.md](SESSIONS.md), the relevant
experiment's agreement/report, and live jobs. Preserve unrelated edits,
immutable scientific artifacts and active executors.

## Scientific analysis autonomy

Jacob explicitly removed the standing per-analysis approval gate, most recently
on 2026-09-18: “remove that AGENTS.md thing bro. I thought we got rid of that.”
His earlier full discretionary permissions for improving this discriminator
remain in force. Declare and record contained experiments, execute them, and
report substantive findings. Do not pause for another approval of each analysis
or routine preparation repair. Stale injected copies of the 2026-09-12 rule and
older pending-approval notes do not reinstate that gate. Ask only for genuinely
missing information or work outside the authorized scope. Preserve production,
scientific integrity, immutable experiments, scheduler rules and concurrent work.

## Earlier water work — 2026-09-18

**Joint occupancy arrangement experiment complete.** Jacob: “I fully approve.”
Read [report](diagnostics/hydration_occupancy_20260918/REPORT.md) and
[commands](diagnostics/hydration_occupancy_20260918/COMMANDS.md). MACE1201908 and
DFT1201910 completed all12patterns/24endpoints (20new,4reused). All16 additions
shift the electronic contrast towardLa; water identity and coupled effects
matter.49fixture tests pass. Occupancy free energies/probabilities remain
unavailable; no new accuracy/default/PQQ claim. Do not rerun completed jobs or
submit source-only native templates. The subsequent local-motion check is complete; see the current result above.

**Water preparation now improves the consumed alpha/GGR discrimination case.**
Read [the result](diagnostics/hydration_network_20260918/REPORT.md) and
[operations](diagnostics/hydration_network_20260918/COMMANDS.md). Exact MACE water
rotations in expanded physical contexts, then original-core/original-recipe DFT,
change alpha−GGR(1GLG) from−14.56/−17.54 to+10.92/+14.57kcal/mol. Reused GGR
replicates give0/6→6/6 directions, one biological comparison, weakest margin0.655.
No threshold fit, altered water count or production/PQQ change.40fixture tests
pass, including archived PQQ energies/bands. Core-transfer job1201867 is complete.

Native orientation1201824 is complete;1201825 remains live. Do not duplicate
or interrupt it. The corrected collector handles Opt's absent final.engrad
and retains actual gradient/geometry pairs; a missing artifact alone is not a
reason to rerun chemistry. Small frozen-coordinate drift remains an explicit
geometry check. Completed MACE jobs1201847/1201849 and DFT adjudication1201853
must not be rerun. The practical route took96GPU-s plus118s on64CPUs for two
structures; the ongoing native searches are separate development cost.

Jacob approved this work with “proceed!”; standing autonomy applies. Collect the
native comparators once they finish. The joint occupancy table is complete; its
separate local-motion follow-on is described above. The old36-task DFT-only manifest
is unsubmitted; do not blindly launch it. Bulk-water reference exists, but missing
bound-state free-energy terms are not zero and no occupancy probabilities exist.
Vault/result email updated. No production promotion or broad validation claimed.

The preceding [hydration-square pilot](diagnostics/hydration_square_20260918/REPORT.md)
is complete:14 endpoints, substantial water-specific effects, one biological
group, no occupancy or accuracy gain established. Do not rerun it or select a
favorable deletion from its outcomes. Recover the current checkpoint for status.

Completed work is preserved in dated reports:

- [Archived-feature classifier](diagnostics/site_classifier_20260918/REPORT.md):
  DFT and DFT+structure retain 25/25 grouped PQQ calls; adding MACE gives 24/25.
  No demonstrated accuracy gain and no out-of-group direct-affinity evaluation.
- [Frozen MACE/DFT comparison](diagnostics/mace_pqq_utility_20260918/REPORT.md):
  MACE retains 25/25 calls at 1.859x median speed; a separate DFT reproduction
  exception prevents the strict combined qualification. Accuracy remains the
  priority; matching consumed references does not establish broader utility.
- [Delegated benchmark augmentation](diagnostics/benchmark_augmentation_20260918/REPORT.md):
  8GY2 prepared as an additional Ca-associated PQQ structural control; no energy
  run or new direct-affinity group. Do not treat structural association as an
  affinity label.

Historical goal and pilot records explain prior experiments; they do not impose
a new approval gate or instruct agents to restart completed campaigns. Recover
current status before execution and avoid rerunning work without scientific value.

Keep the existing baseline/default accessible. Buffered inbox PQQ, canonical
fixed-core PQQ, generic v2, repaired peptide v3 and environmental challengers
are distinct protocols. Use only a protocol's own reference/decision policy;
missing or failed calculations remain explicit, with no baseline substitution.

Follow [DIRECTORY_POLICY.md](DIRECTORY_POLICY.md): candidate compute products
go under `workspaces/`, compact experiment records under `diagnostics/`.
Use the existing pinned runners and environments. Commit only your own scoped
changes and record the completed work and remaining live work in SESSIONS.md.
