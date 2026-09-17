# Active MACE discriminator work — 2026-09-17

**Goal active.** Jacob authorized autonomous analyses and available resources;
the old per-analysis permission gate has been removed from home/project AGENTS.
Preserve the production baseline, existing results and other agents' jobs.

## Completed whole-chain candidate: rejected for promotion

`mace_omol_intact_chain_matched_coordination_v1`: pinned OMOL100M, whole prepared
chain, fixed geometry/microstate/waters, native vacuum energies, no forces.
R_coord=(E_bound,Ca−E_detached,Ca)−(E_bound,La−E_detached,La). Larger values
are relatively more La-like; zero and old ORCA bands have no authority here.

The five-case development panel passes all three declared relative comparisons:
XoxF−MxaF+9.03kcal/mol; alpha−GGR+37.15/+28.09 for two alpha structures.
One qualified affinity comparison and a separate PQQ functional-class test;
all consumed. It subsequently failed the full canonical calibration and the
disconnected-spectator consistency test. No default promotion. See
[result](../mace_omol_20260917/INTACT_REPORT.md).

All25canonical cases computed, with a negative separation gap−306.006296kcal/mol.
No classification bands were issued; two transfer raw scores are present and
1KB0 is unsupported. See [canonical report](../mace_omol_20260917/INTACT_CANONICAL_REPORT.md).
The sodium test shifts four scores4–15kcal/mol and reverses XoxF/MxaF. Exact
checkpoint inspection finds identical features for all charges−100..−6,
explaining MxaF's apparent invariance. See [spectator report](../mace_omol_20260917/INTACT_SPECTATOR_REPORT.md)
and [embedding audit](../mace_omol_20260917/CHARGE_EMBEDDING_AUDIT_REPORT.md).

Exact edge/product batching is qualified against native CPU energies. The
whole-chain model fits one A5000; recent PQQ endpoints take roughly36 seconds
including model loading. This is a local learned model with global charge
conditioning, not a demonstrated full electrostatic description. The
[locality audit](../mace_omol_20260917/INTACT_LOCALITY_REPORT.md) establishes
cancellation of distant and additive embedding terms.

## Execution state

| Job | Purpose | Inventory |
|---|---|---|
|1200819|Completed; calibration failed|104new endpoints+8exact reuses|
|1200823|Completed; consistency failed|20new endpoints+20source references|
|1200828|Charge-feature ablation submitted|42new descriptor forwards|
|1200809|Original native H200 numerical reference|Earlier unchanged14-task manifest; pending|

Do not duplicate or alter these jobs/manifests. Check live state before acting.
Their workspaces are under `workspaces/mace_omol_20260917/`:
`intact_panel_run_v1`, `intact_spectator_v1`, `intact_qualification_v1`.

**1KB0 whole-chain preparation is invalid:** two false peptide connections
cross missing structure. Its raw endpoints are diagnostic-only. The mandatory
audited canonical reporter is frozen under
`intact_panel_reporting_source_v2/implementation/`. Follow the exact
[commands](../mace_omol_20260917/INTACT_COMMANDS.md), retaining1KB0 as unavailable
in the three-transfer denominator. Old fixed-core baseline1KB0 is unchanged.
New strict preparation v2 gives27supported/1unsupported; original run inputs
remain immutable. All five earlier development preparations separately pass
the same1718-bond audit (`intact_five_integrity_v1/result.json`).

A single bounded local waiter completed:
`workspaces/mace_omol_20260917/intact_panel_report_after_1200819.py`.
It invoked the exact audited v2 reporter once, with timing receipts. The report
is complete under `intact_panel_report_v1/`; no waiter remains active.

The spectator construction adds one explicitly constructed distant Na+ while
retaining every original coordinate and protein formal charge. Native MACE
only sees the total charge; it cannot certify that fragment charge assignment.
Read its [frozen scope](../mace_omol_20260917/INTACT_SPECTATOR_PLAN.md). The
existing runner collects numerical and consistency gates separately. A failed
consistency gate is not a failed execution or a biological misclassification.

## Active next model

[CHARGE_ABLATION_PLAN.md](../mace_omol_20260917/CHARGE_ABLATION_PLAN.md) declares
one fixed representation change: mask the raw global charge embedding with
zeros before the native joint projection, preserving spin, weights, physical
charges and every coordinate. Its outputs are explicitly energy-like learned
descriptors, not quantum energies for the recorded electronic state.
Protocol `mace_omol_intact_charge_feature_ablation_descriptor_v1`.

Initial42forwards:8core native/batched qualification,14alpha numerical checks,
16other whole-chain primaries,4GGR sodium checks. Only a passing numerical and
three-direction gate permits100canonical forwards with8qualified crystal reuses.
No molecular ablation endpoint has run yet. The adapter is implemented and a
CPU component check passed under `charge_ablation_component_v3/`:201charge
categories give identical joint features, matching the explicit zero-charge plus
original-spin projection exactly; no parameters changed and no molecular call ran.
Cost17.27wall/21.53CPU seconds,1616720KiB peakRSS. V1failed from an import-name
collision; V2test incorrectly assumed feature concatenation order. V3uses the
checkpoint's actual spin-then-charge order. All attempts remain; no model setting changed.

Component qualification and runner integration are complete. The exact42-task
manifest is submitted as1200828 under `charge_ablation_development_v2/`, SHA
78c781671295658656111f94f9b7cd6625b5bb696ab85e726fe70ea5de2d30c5.
Four preflight tests and six native regressions passed. V1never ran; v2fixes
receipt validation for the actual atom-expanded charge-feature row count.
The runner requires core qualification before whole-chain calls and alpha
numerical qualification before the remaining primary/spectator tasks. Do not compare these
outputs to native OMOL as if the model were unchanged. Earlier physics failures
must remain visible even if this new descriptor predicts well.

## Evidence and remaining work

Both reports retain actual receipts, denominators and measured costs.
Do not change criteria after inspection. Preserve numerical credibility,
predictive usefulness and affordability as separate judgments. Earlier failed
POLAR, mechanical-response and core descriptor experiments remain recorded;
neither their failure nor this five-case success completes the goal.

Canonical calibration tests breadth within PQQ. It cannot prove added value
beyond composition or general affinity discrimination. The evidence gap for
new independent non-PQQ directions remains in the existing challenge curation;
do not turn protein-level measurements into independent site labels.

Own implementation commits through canonical preparation:19cede4,ed2e65c.
Later locality/spectator changes and the declared ablation are scoped separately. Many other
repository edits belong to concurrent/historical work; never blanket-stage.
