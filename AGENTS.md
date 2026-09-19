# Agent entry point

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

Jacob approved water-network/occupancy modeling with “proceed!”. Recover
[the hydration-network plan](diagnostics/hydration_network_20260918/AGREEMENT.md).
Native DFT water-orientation searches are running as jobs 1201824 (1F6S) and
1201825 (6IP9), four searches per job on64CPUs, two seeds for each metal.
Only inner-water H orientations change; protein and water oxygens stay fixed.
New70/76-atom cores include source-defined missing hydrogen-bond neighbors;
protein composition is common across the two structural replicates. A separate
gas-water reference/comparator completed as1201831 after a launcher-only
retry; see WATER_REFERENCE.md for results and actual whole-node allocation cost. No occupancy probabilities,
new classification or production change yet. Preserve current jobs and collect
actual outcomes before extending the finite occupancy-state manifests.

The native MACE proposal check passed (1201847); all eight exact rigid-water
MACE optimizations subsequently converged (1201849). Four selected configurations
are undergoing native DFT adjudication as1201853. Read MACE_PROPOSAL_REPORT.md and
MACE_OPTIMIZATION_PLAN.md. Preserve the native DFT comparison jobs. No broader
accuracy/occupancy gain is established yet. The corrected Opt collector retains
actual printed gradients and explicitly handles the missing final.engrad;
do not rerun these searches merely for that artifact. Frozen-coordinate drift
in native Opt remains an explicit geometry check, not an erased failure.

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
