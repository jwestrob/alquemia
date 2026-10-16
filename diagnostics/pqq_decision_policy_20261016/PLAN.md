# Calibration-only PQQ decision experiment

Declared 16 October 2026 before new threshold evaluation. Jacob approved keeping
the chemistry fixed and testing whether the decision rule withholds useful calls.
This supersedes the discussion pause for this saved-data comparison only. It does
not resume staged-executor development or authorize production promotion.

## Question

Do thresholds learned solely from the designated 25 canonical proteins resolve
current abstentions without introducing wrong calls on the consumed structural
benchmark? Global class ordering has already been inspected, so this is explicitly
retrospective development evidence. Evaluation extrema must not set thresholds.

## Fixed data and methods

- Use the completed corrected-parameter, frozen-charge common6 scores unchanged.
- Use the original canonical25 IDs from the frozen-distribution inventory, their
  common6 canonical collection, and all257 unique current transfer contexts.
  Keep260 memberships, three aliases and18 superseded records distinguishable.
- Primary analysis uses operational selection; repeat the same rules for the
  mathematical-selection scores. Do not select a favorable policy afterward.
- Use supplied expected_class labels, never conditioning metal as a label.
- Deduplicate by context_id and define the noncanonical subset by excluding every
  context belonging to a designated canonical record, including aliases.
- Retain all100 existing triple definitions, median-score aggregation, and six
  preparation exclusions. Do not replace median aggregation with voting.

## Rules frozen before evaluation

For each selection policy, using canonical scores ONLY, let C be maximum Ca score
and L minimum La score. Require both classes, complete calibration and C < L.

1. **Existing bands:** exact archived Ca_max/La_min decision rule, reproduced as
   the comparator. No change to its inclusive boundaries or recorded outcomes.
2. **Primary midpoint:** T = C + (L-C)/2. Predict Ca for score < T and La for
   score > T. An exact tie remains inconclusive. No margin is fitted.
3. **Secondary calibration-sensitivity band:** remove each biological_group in
   turn from the canonical set and obtain its midpoint T_minus_group. Set the
   lower/upper bounds to min/max over T and all those group-omitted midpoints.
   Predict Ca strictly below the lower bound and La strictly above the upper
   bound; otherwise inconclusive. No multiplier or width search. This is an
   empirical leverage diagnostic, NOT a confidence or prediction interval.

If a required fit lacks either class or fails separation, that rule is unavailable;
no transfer-informed fallback. Missing scores remain unavailable. Report every
rule, regardless of outcome; the midpoint remains the primary test.

Separately report leave-one-biological-group-out predictions on the held-out
canonical group using ONLY its T_minus_group. Do not call full-envelope predictions
held-out: other envelope constituents include that group's data. This diagnostic
does not make historically consumed model-development cases a blind test.

## Outputs and checks

Report correct/wrong/inconclusive/unavailable outcomes for all257, canonical and
noncanonical subsets, source-conditioning strata and original cohorts. List every
changed decision, especially resolved abstentions, newly wrong and lost-correct
calls. Report per-protein outcomes and the original100 triple medians. No new
source-selection or ensemble rule is introduced. Report reference leverage and
the original abstention cases under every rule without retuning.

Verify exact archived-band replay, canonical membership/score/label joins,
deduplication, unchanged medians and exclusions. Algebra checks use actual pinned
calibration records: midpoint within gap, exact-tie abstention, held-out removal,
and refusal of explicitly incomplete copies. No fabricated molecular data.

## Execution and ownership

One stdlib-only job on one shared Slurm CPU, no memory/time override and zero
MACE/GFN2/DFT calls. Root owns source, submission, report and integration. Separate
agents inventory the designated records and review the frozen rules/source.
Reuse existing completion watcher. Do not act on stale campaign notifications.
Compact code/plan/results live here; execution products go to the matching
workspaces directory. Freeze input/source hashes before submission. No expensive
execution framework, plotting package, external survey or new molecular campaign
is needed to answer this question.
