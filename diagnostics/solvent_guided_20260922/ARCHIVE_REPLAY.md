# Archived geometry rankings justify a small solvent-guided search test

The current pipeline already selects its final geometry using the composite
MACE+GFN2(ALPB−vacuum) score. The archived replay shows why optimizing proposals
only under vacuum MACE can miss the geometry preferred by that final score.
It does not establish an additional accuracy gain from the proposed new search.

| Sources | Complete / all | Ca choice disagreements | La choice disagreements | Ca / La regrets >0.1 | Sources with absolute ΔR >0.1 |
|---|---:|---:|---:|---:|---:|
| Original canonical/crystal/unknown30 | 30/30 | 14 | 6 | 13/4 | 17 |
| All noncanonical folds | 204/225 | 93 | 29 | 90/24 | 112 |
| Ca-conditioned folds | 107/125 | 47 | 15 | 45/12 | 55 |
| La-conditioned folds | 97/100 | 46 | 14 | 45/12 | 57 |

Both choices use the same available common geometries and the existing0.1model
kcal/mol origin-retention policy. Regret is the composite energy at the
native-selected geometry minus the best composite energy in that same pool.
Across complete folds, maximum regret is32.35656kcal/mol for Ca and4.82686 for La.
The Ca outlier is Q92WY9 La-conditioned sample0; it remains in the summary and
was not added to the predeclared common8 pilot to chase its large effect.

Composite-selecting rather than native-selecting geometries changes R by
−32.35656 to+3.82538kcal/mol across the204folds. The median iszero: many sources
already agree. This is a difference between selections **within the existing
computed candidate pool**, not a new biological score or a result of new calls.
Mathematical minima are retained separately from operational origin retention.
No classification threshold or label entered ranking.

All21unavailable fold sources remain:20inherited base-pool failures and one
unavailable La adaptive proposal. Source groups use the recorded root-case to
biological-group mapping; canonical/fold instances of the same protein are
merged. Crystal controls and unknown PLM sources retain separate record roles.
The225folds are repeats of25consumed reference proteins, not independent biology.

The exact source/metal choices, candidate work components, regrets, signed ΔR,
conditioning summaries and per-group summaries are in
`workspaces/solvent_guided_20260922/archive_rankings_v2.json`.
Version1 remains intact; version2 only unifies canonical/fold biological-group
identifiers through the recorded root-case mapping. Numeric results are unchanged.
No molecular calculation, training, parameter fit or threshold change ran here.

Next: execute only the frozen [finite common8 probe plan](PROBE_PLAN.md), after
parent integration review of the actual geometry manifest. This probes for new
useful candidates; no full225 scoring is automatic.

## Subsequent numerical caution

During the parallel source-start track, a near-identical Q88JH5 geometry produced
a several-kcal/mol native GFN2 solvent-contrast discrepancy while MACE agreed.
That audit is separate and ongoing. Archived regrets above remain literal
computed-energy differences, but cannot yet be uniquely attributed to physical
solvent response or used to claim predictive improvement. The finite point plan
and already-running calculation set remain unchanged; no numerical rescue or
new threshold is introduced here.
