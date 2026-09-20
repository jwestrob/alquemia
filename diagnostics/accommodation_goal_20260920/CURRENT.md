# Overnight accuracy research checkpoint — 2026-09-20

## Latest substantive result and continuation

Native DFT confirms both large PLM **differential** donor responses, not just
their calcium components. MACE captures the direction and approximate magnitude;
the solvent term worsens differential agreement. All 16 original native endpoints
completed; all 12 planned work comparisons are available. This supports a physical
accommodation mechanism, not known PLM labels or validated relaxed scores.
See [final native report](../accommodation_torsion_20260920/FINAL_DFT_REPORT.md).
Existing profiles have boundary minima/nonconstant curvature. The original native
job1203771 is complete; do not rerun it. Final native result/figure are saved.

The [next declared experiment](../accommodation_nonlinear_20260920/PLAN.md) tests
actual composite energies/analytic gradients along one or two whole-carboxylate
torsions. Second_shell's implementation4755719 passed preflight and runs as1204162.
The first Ca control is an interior qualified cheap minimum; the corresponding La
candidate is stationary but fails the frozen curvature refinement tolerance.
That is a numerical qualification failure, not evidence of an unstable basin.
All 28 exact reference mappings are prepared (e6ef185), with no energy calls.
Root owns technical review, execution coordination and later native validation.
The first experiment contains eight metal-specific starts on the four original
contexts. No fresh DFT is included until a separate validation manifest exists.
No score-directed starting point, fitted spring or harmonic entropy is introduced.
Water_basins implemented the separate native candidate adapter82dd519; it still
requires final pilot output and genuine task-path preflight before submission.
Root identified a shard-directory containment fix before any candidate calculation.

Khoury owns the new [proposal experiment](../accommodation_nonlinear_20260920/PROPOSAL_PLAN.md):
all 28 reference contexts plus the two original PLM contexts. Native MACE proposes
the bounded donor arrangement; the existing primary composite evaluates the
proposal and original, selecting only by energy. This tests whether solvent must
be solved throughout the search. It preserves the running full-gradient pilot,
does not remove solvent from scoring and does not promise a continuous minimum.
Its source compatibility/preparation is in progress; no submission yet at this
checkpoint. No threshold fitting or production/default change.

The [latest partial benchmark](../accommodation_fold_DFT_20260920/PARTIAL_MATCHED_REPORT_v2.md)
contains189/416 completed DFT endpoints. On the same42 La-conditioned source pairs,
context composite corrects one additional call and removes the one DFT wrong call;
on37 available triples DFT, native core and context composite all agree correctly.
This is a small consumed-structure robustness result; completion order and missing
coverage prevent a panel conclusion. The four baseline DFT jobs continue normally.
The first PLM response observer completed its read-only task; it launches nothing.

Third substantive email accepted by local relay at12:10UTC; exact message/receipt
under workspaces/accommodation_goal_20260920/email_physical_response_v1/. Earlier
dated checkpoints below are retained as history, not new execution instructions.

## Why these experiments

The release now makes source-backed PQQ scoring practical. The active goal goes
further: stop a single imperfect model from controlling a metal-specificity call,
and identify a physical adjustment that improves model reliability. The saved
reference folds test structural robustness against known classes. Whole-donor
rotations test whether diagnosed compression can be relieved with a physically
credible metal-dependent response. Unknown PLM predictions remain unlabeled.

## Current execution

**Completed benefit:** on the additional La-conditioned reference folds, the
context solvent correction removes the native MACE wrong calls. The predeclared
four-fold median also corrects Q9Z4J7 on unchanged protein coverage. This supports
more reliable use of La-conditioned folds; it does not establish superiority to
DFT or new biological accuracy. Ca-conditioned folds expose two new individual
errors, and incomplete preparations still limit coverage. Read the
[full comparison](FOLD_COMPARISON_REPORT.md), committed as7069cd1.

The next practical question is whether this robustness survives the three-fold
inventory available for the PLM scan. The separate [three-fold plan](THREE_FOLD_PLAN.md)
tests every possible three-of-four subset of the consumed La-conditioned reference
folds, with unchanged bands and no new molecular calculations or selected subset.
That addendum is complete (3048586): all available context-composite triples are
correct, and native core matches that result on identical coverage. Read
[three-fold report](THREE_FOLD_REPORT.md). It supports subset robustness, not
superiority to every simpler method or validation of PLM labels.

Second_shell now owns an optional `standard ensemble` aggregation interface:
three declared La-conditioned source folds, strict membership/state/sequence
checks, individual scores/spread/median and developmental frozen-band transfer.
One real end-to-end request will use the lexicographically first three
noncanonical Q9Z4J7 La source IDs; six MACE and twelve GFN2 evaluations, no DFT.
This validates the interface on consumed sources, not new predictive evidence.
The interface check makes no ensemble default or production scan changes.

The reference interface check completed successfully as1204055: all declared
source preparations/endpoints succeeded and the three-fold aggregate reproduced
archived values. The optional report retains individual scores, spread and
unknown-input-domain flags. Its release/ensemble checks pass (23tests), committed
as11d9669. Second_shell owns a subsequently authorized isolated continuation on the TWO existing PLM compression
diagnostics, all three saved folds each (12MACE/24GFN2, no DFT/folding). This asks
whether their predictions depend on one source fold; it does not modify the
production scan or supply biological labels. See [source readiness](PLM_SOURCE_READINESS.md).
One previously inadmissible source stays included and flagged; preparation
failure must produce an unavailable group, never a replacement member.
Its separate plan is diagnostics/plm_fold_sensitivity_20260920/PLAN.md;
independent jobs1204068(07ab) and1204069(8344) use the existing H200 layout.
Both are now complete: all six preparations and12MACE/24GFN2 endpoints succeeded.
Every fold of both targets remains Ca-supported; each three-fold median equals
sample0. Fold aggregation therefore does not resolve this particular concern.
Extra-Asp compression persists across all six structures. Sample0 contexts,
states and methods match the original torsion origins, with exact/numerical replay.
The old admission failure is retained even though current source preparation
supports that member. Biological labels remain unknown; do not call this six
classification failures or alter the reference bands to make them La-like.
The native donor-response validation is the next relevant evidence for this
persistent geometry concern. No extra angles or replacement samples were added.

Khoury owns the same100-triple comparison against preserved DFT. Its partial
snapshot is insufficient for a panel conclusion; use the ready final command
only after the existing DFT jobs finish. No additional DFT work is needed.
Comparator/tests/report are committed as11a41b0, with five actual-artifact tests
passing. The ready command is in accommodation_fold_DFT_20260920/THREE_FOLD_DFT_REPORT.md.

- Fast-PQQ release complete: commit8dd81ba, report in pqq_fast_release_20260920.
- Matched source inventory complete: commit8ace781. All250 existingAF3 structures
  pinned, including125Ca-conditioned and125La-conditioned, no new folds.
- Root source preparation1203744 complete:233supported,17unsupported after retaining
  seven machine-last-bit replay differences as numerically identical. Only the
  comparisons are reconciled; no coordinates changed. Original failures persist.
- All834fresh native MACE endpoints complete in76.214s in a warm calculator;
  98canonical endpoints reused. The2fresh execution verification endpoints exactly
  match their archive. No scientific MACE failure; earlier1203747 failed importing
  CPU-only gemmi before any scientific evaluation and was repaired with lazy imports.
- Solvent tasks:1668new +196exact-compatible reused nativeGFN2 endpoints.
  Four disjoint case manifests completed as1203797/1203798/1203799/1203800 on existing64CPU runners.
  Two vacuum SCF failures remain unavailable; all other results were collected.
  Initial split failed path-containment validation before execution; version2
  stages actual input copies within each manifest directory. No solver results
  were generated or overwritten in either preparation failure.
- Separate physical profiles:64MACE endpoints and126/128successfulGFN2 outputs.
  Relieving compressed extraAsp donors preferentially lowersLa energy in both
  diagnosed PLM examples; this is not biological accuracy evidence. Most of the
  response is in OMOL, with the solvent term moderating it. NativeDFT1203771 now
  checks the16predeclared endpoints. Keep the failed1H4IGFN2 point unavailable.
  The first four DFT control endpoints now agree with all displacement-work signs
  and small differential response; the PLM differential response remains pending. See
  [partial DFT report](../accommodation_torsion_20260920/PARTIAL_DFT_REPORT_v1.md).
  Partial_v2 now adds the first PLM Ca origin and both displaced Ca endpoints:
  native DFT confirms the large predicted energy penalty/relief in both directions.
  This supports actual local strain response on the calcium surface. La response
  and metal selectivity remain unavailable; do not infer them from calcium alone.
- CC1202429 remains independent; preserve its monitor and outputs.

## Numerical replay and pooling rules

Source/PQQ/charge/water selection stays fixed. The seven canonical context
comparisons differed only by2.22e-16Å in constructed cap coordinates on another
CPU, with identical normalizedPDB bytes. The reconciliation accepts <=1e-12Å,
far below source3-decimal/core6-decimal precision, before inspecting scores.
It records raw comparisons and original statuses and performs no coordinate edit.
All other preparation failures remain unsupported.

Original25canonical geometries are replay controls. Primary pools are the other
fourLa-conditioned and fiveCa-conditioned samples per protein. The declared
medians/equal-weight mean require all members; missing preparations do not become
favorable partial ensembles. Core/composition invariants and changing context
membership remain separate. Frozen old bands are a transfer test, never refitted.

## Current files

Products: workspaces/accommodation_goal_20260920/folds_v1/.
Sourcepreparation/preparation.json is raw; preparation_reconciled_v1.json retains
last-bit acceptance; mace_v3/manifest.json is actual successful execution;
INVENTORY.json is its233-case native endpoint inventory. mace/ and mace_v2/
contain unexecuted or import-failed preparation attempts, not successful science.
The unexecuted full solvent manifest is partitioned only for throughput; use
solvent_shards_v2/ after validation and submission, not the full parent manifest.

Three actual-source adapter/reconciliation tests pass. Additional report tests
are owned by second_shell. No production rescore, PLM output change, new threshold,
mechanical correction, or ensemble promotion has occurred.

## Contained continuations

Geometry-only mapping41851ec found62flagged reference folds, but only one has
the same extraAsp compression as the two PLM examples. The others are anchorGlu
warnings and do not inherit the extraAsp validation. The sole matched-phenotype
reference (MMOL1770,Ca-conditioned sample1) completed fixed±0.2rad with4newMACE
and8newGFN2 calls, q0reused, in accommodation_reference_profile_20260920. No newDFT.
Both metals prefer the same boundary point: native MACE's wrong call becomes
inconclusive; the already-correct composite call gains margin. This does not add
a correct decisive classification or establish an interior relaxed minimum.

A separate baseline comparison is running as1203976–1203979: exact old
DFTcore recipe on all208supported noncanonical fold samples (416endpoints), with
25canonical outputs reused and17unsupported retained. This is a direct comparator
for accuracy/conditioning robustness, parallel to candidate development, with no
threshold fitting or baseline change. Four64CPU jobs each run four16-rank endpoints;
automatic collection is installed. Read
[status and commands](../accommodation_fold_DFT_20260920/STATUS.md), commit82c7bcc.
Do not duplicate these jobs. Early partial pairs are insufficient for a panel claim.

The H-only minimizer diagnosticb62b011 converged under the same generic potential
after monotonic preconditioning. It fixes silent numerical failure; that potential
still produces longHbonds, also present in its nominally successful controls.
No Hmodel change is promoted. See accommodation_controls HYDROGEN_DIAGNOSTIC_REPORT.
