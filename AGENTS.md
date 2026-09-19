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

The matched observed-water mechanism test is complete: job 1201801, all 14
ORCA endpoints successful. Repairing water H geometry shifts the alpha contrasts
+4.93/+2.48 kcal/mol toward La. Individual retained waters shift them from −19.95
to +5.01 kcal/mol; hydration identity matters, but occupancy and improved accuracy
are not established. See [report](diagnostics/hydration_square_20260918/REPORT.md)
and [checkpoint](diagnostics/mace_discriminator_goal_20260916/CURRENT.md).
No hydration job remains live. Baseline and old inputs remain unchanged.
Next research direction: explicit, consistently referenced hydration states;
do not select deletions by their agreement with labels or rerun this pilot.

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
