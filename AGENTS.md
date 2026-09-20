# Agent entry point

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
