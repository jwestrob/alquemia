# Agent entry point

## Current: local water response passes — 2026-09-18

All four metal states (1F6S11 and6IP9110, Ca/La) pass the prescribed local radial
energy/gradient/curvature checks. Eight native DFT evaluations and16MACE calls
completed as1201953/1201958; startup-only1201954 failed before inference and is
included in cost. Curvature errors0.58–2.53%; maximum anchored energy error
0.00707kcal/mol and Ca−La response error0.00905kcal/mol. This supports continuing
DFT-anchored cheap water mechanics; it does not yet supply occupancy/entropy or
new classification accuracy. Baseline/default/PQQ unchanged.56real-fixture tests
pass, zero skips. Cost:69,008allocated core-s and137GPU-s. No next-stage run.

Read [motion report](diagnostics/hydration_motion_20260918/REPORT.md) and [commands](diagnostics/hydration_motion_20260918/COMMANDS.md).
Do not rerun completed motion or occupancy jobs. Final result RESULT_v2.json and
export_v3 use a frozen analyzer. Native orientation comparator1201824 has completed all four1F6S optimizations.
They converge normally but fail the frozen-coordinate tolerance (drift2.26e-5 to
7.46e-4Å); one also exceeds the rigid-water tolerance. Actual final energies and
gradient traces are retained in orientation_1f6s_v1/collected_opt_v2.json, without
promoting them to qualified exact-geometry minima. No rerun or tolerance change.
Comparator1201825 for6IP9 remains live; preserve it. These are separate older
runs and do not affect the successful exact-coordinate motion checks.
Next scientific target: coupled physical water motions and stable basins.


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

## Current work — 2026-09-18

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

Jobs1201824/1201825 remain live native DFT orientation comparators; do not duplicate
or interrupt them. The corrected collector handles Opt's absent final.engrad
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
