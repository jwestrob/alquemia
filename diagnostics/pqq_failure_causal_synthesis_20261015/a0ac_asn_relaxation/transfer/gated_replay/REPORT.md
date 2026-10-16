# Screening retains the full amide-search result on the consumed archive

The retrospective gate reproduces all257 full-search scores, both selected metal
geometries, decisions and availability exactly, under both reporting policies.
Only A0ACD6B9F2 Ca-conditioned sample4 and La-conditioned sample4 trigger.
This supports implementing a live opt-in candidate; it is not a live speed test
or evidence of recall on unseen geometries.

## Fixed rule

Use the original common4 operational composite selection, including its strict
0.1kcal/mol improvement requirement and original-first tie handling. If either
metal selects `amide_flip`, admit both bounded relaxations and the common6
cross-scored pool. Otherwise retain common3. An unavailable screen or required
pool yields unavailable. Both operational and mathematical reporting follow the
same operational trigger. Routing reads availability and candidate identities;
class, source conditioning and desired sign do not select the branch.

The rule was declared after inspecting the archive. It is label-free at execution,
but retrospective in design. The short-contact fixed flip is a screening geometry;
this result does not independently validate it as a relaxed physical endpoint.

## Actual results

| Quantity | Result |
|---|---:|
| Unique contexts / cohort memberships | 257 /260 |
| Triggered / untriggered / unavailable screens | 2 /255 /0 |
| Correct / wrong / inconclusive / unavailable | 249 /0 /8 /0 |
| Score, selection, decision or availability differences from full common6 | 0 |
| Eligible paired-search slots / avoidable slots in replay | 4 /510 |
| New molecular calls | 0 |

La selects the fixed flip in both triggered screens; Ca selects `adaptive_Ca`.
The one unique classification repair remains the Ca-conditioned sample4, shared
by both cohorts. All260 memberships, three alias contexts and18 superseded
historical records remain accounted for. Canonical bands are unchanged.
Independent post-reboot review checked20 inspected source/result pins and found
no discrepancy. [Compact machine-readable evidence](RESULT.json) pins the full
workspace result, including all individual routed records.

## What this does and does not buy

On these saved inputs, screening preserves the useful relaxed geometry and
reduced source spread without requiring every context to enter the extra search.
The510 avoided slots are policy eligibility counts, including historically reused
contexts. They are not completed jobs saved, recovered compute, or measured time.
The screen itself costs molecular evaluations. A fixed alternative that initially
loses can relax into a useful basin; this archive cannot establish the screen's
recall elsewhere. The eight remaining abstentions and electronic/reference
limitations in the [full report](../FINAL_REPORT.md) remain.

Next: a separately declared live candidate scout using existing executors, with
both triggered contexts and an untriggered canonical control, before claiming
runtime benefit. No new molecular scout or production promotion occurred here.

## Execution and reproduction

Job6045 completed with exit0:0 in2 scheduler seconds on1CPU; application elapsed
0.9101868390571326s. The allocation was normal shared Slurm, with zero GPU time
and zero new scores, optimizations, MACE, SCC, native solvent or DFT calls.
The unchanged saved-summary replay exactly matched its previous result before
routing. Full result SHA256:
`9e81ca26276015aaf51ee24bed58f5b1e8350725823b761a399482989c3c3425`.

[READY.md](READY.md) retains the exact executed command and source pins;
[PLAN.md](PLAN.md) is the original immutable rule declaration.
Workspace: `workspaces/pqq_failure_causal_synthesis_20261015/asn_transfer_gated_replay_v1/`.
The reboot did not invalidate these completed artifacts. No rerun was performed.
