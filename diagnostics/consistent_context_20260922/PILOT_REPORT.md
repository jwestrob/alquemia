# Fixed-context calibration and initial fold transfer

Status: canonical 25 / crystal 3 and pilot 36 complete; unchanged full 225 continuation
is being prepared. These are consumed development sources, not new biological
validation. The released scorer and its reference remain unchanged.

## Initial outcome

The new canonical-only gap is **11.066622 model kcal/mol**, compared with the
released static-context gap of 5.076366. All 25 canonical and all 3 crystal controls
classify correctly. Crystals are exact singleton identities with their old
component receipts reused. The reference was frozen before transfer execution.

On the four-group pilot, static context gives 21 correct, 2 wrong, 2 inconclusive;
fixed membership gives 19 correct, 1 wrong, 5 inconclusive. The denominator is 36:
11 previously unsupported structures remain unavailable in both methods. Thus
one real error is corrected, while three previously correct calls become
inconclusive. This is **no demonstrated net accuracy improvement**.

All strict La4/Ca5/equal-arm decisions stay unchanged. Across all 16 declared La
triples,13 correct / 3 unavailable becomes11 correct / 2 inconclusive / 3 unavailable. The
same simple native-core method already gives13 correct / 3 unavailable.

## Changed single-source decisions

| Source | Static context | Fixed membership | Raw R shift, model kcal/mol |
|---|---|---|---:|
| a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-0 | correct | inconclusive | -1.004855 |
| a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1 | wrong | correct | 8.740205 |
| a0a3f2yly8-pqq-la_model__conditioned_La__seed-1_sample-0 | correct | inconclusive | -1.309832 |
| a0a3f2yly8-pqq-la_model__conditioned_La__seed-1_sample-3 | correct | inconclusive | -1.037021 |

The corrected A0A3F2YLY8 Ca-conditioned sample 1 shifts +8.740205 overall, composed
of −5.731701 native-OMOL and +14.471905 solvent-transfer contributions. The
expanded preparation changes composition/cavity along with retained surrounding
chemistry; this does not isolate a mechanical or hydrogen-bond cause. The three
new inconclusives shift about −1 kcal and encounter the independently calibrated
wider inconclusive interval. No band was moved after those outcomes.

Within La-conditioned samples, score ranges change 21.56→28.89 for A0A3F2YLY8,
22.75→19.69 for A0ACD6B9F2,21.91→15.97 for Q9Z4J7 and6.58→9.74 for Q88JH5.
Fixed composition does not uniformly reduce structural sensitivity.

## Executed calls and preserved records

Canonical stage:30 fresh MACE +60 fresh GFN2;26 MACE +52 GFN2 exact compatible
components reused. Jobs 1209972 and1209977 completed. Job 1209973 was cancelled
while still pending with zero elapsed/scientific work because its standard
partition placement was unavailable; the same manifest ran on a CPU-only GPU
partition allocation.

Pilot stage:40 fresh MACE +80 fresh GFN2;10 MACE +20 GFN2 components reused. Jobs
1210013/1210014 completed. All fresh molecular calls in these stages succeeded.
Canonical worker wall 8.541 s; pilot worker wall 9.066 s. These are warm worker
component costs, not source-to-score end-to-end timing. Full allocation records
are in SACCT_INITIAL.txt.

Eleven real-fixture tests pass (TESTS_v2.txt), including original-coordinate and
paired-state invariants, exact old selection/crystal replay, source graph union,
actual component algebra, frozen calibration, strict missing-member medians and
corrupted-reference rejection. No molecular calculation ran as a test.

The first calibration collection/reference reported unavailable because a parser
status field accidentally replaced the collector's completion status. Original
records remain preserved; collection_final_v2 and REFERENCE_v2 fix the status
merge using the same actual energies, with a regression test. No scientific
calculation was repeated for this issue.

See FULL_TRANSFER_PLAN.md for the fixed continuation. All raw inputs, receipts,
components and both comparison tables remain under the separate
workspaces/consistent_context_20260922/ tree.
