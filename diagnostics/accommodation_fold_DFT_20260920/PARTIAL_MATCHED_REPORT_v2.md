# Partial matched DFT/MACE comparison

2026-09-20. One fresh immutable snapshot, requested by root; existing policies,
bands and sources unchanged. No new molecular calls or implementation changes.

**The completed three-fold DFT summaries match the strongest MACE arms.**
Individual folds show real disagreements, but this partial result does not
establish overall superiority. The snapshot has 189/416 DFT endpoints complete,
227 pending and no recorded failures. It supplies 91/225 primary source pairs,
37/100 triples, nine La4 pools, nine Ca5 pools and eight balanced pools.

## Same completed sources

Counts below are restricted to the **same 91 primary sources** with complete
DFT pairs. Another 134/225 primary sources remain unavailable to DFT in this
snapshot; all are retained, including the 17 unsupported preparations.

| Method | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Preserved DFT | 84 | 3 | 4 | 0 |
| Native MACE core | 84 | 0 | 7 | 0 |
| Native MACE context | 87 | 3 | 1 | 0 |
| Composite core | 78 | 0 | 13 | 0 |
| Composite context | 86 | 2 | 2 | 1 |

The three DFT wrong calls are A0ACD6B9F2 Ca-conditioned sample 4 and
La-conditioned sample 4, plus MMOL1770 Ca-conditioned sample 1. All are known
La-family references called Ca-supported by DFT. Native core makes all three
inconclusive. Native context corrects both A0ACD6B9F2 cases but retains the
MMOL1770 error; context composite corrects MMOL1770, retains the
Ca-conditioned A0ACD6B9F2 error and makes its La-conditioned case inconclusive.
Context methods also have different errors among other sources. This is a
tradeoff, not a uniform correction.

### Conditioning-specific results and every DFT noncorrect case

The main scanner uses La-conditioned folds. On exactly the **42 completed
La-conditioned sources**, context composite gives one additional correct call
and removes the one wrong call relative to DFT. This is a small observed gain
on a partial, consumed set; three-fold summaries still show parity.

Counts are correct / wrong / inconclusive / unavailable on the matched rows:

| Method | Same 42 La-conditioned sources | Same 49 Ca-conditioned sources |
|---|---:|---:|
| Preserved DFT | 40 / 1 / 1 / 0 | 44 / 2 / 3 / 0 |
| Native MACE core | 40 / 0 / 2 / 0 | 44 / 0 / 5 / 0 |
| Native MACE context | 40 / 1 / 1 / 0 | 47 / 2 / 0 / 0 |
| Composite core | 37 / 0 / 5 / 0 | 41 / 0 / 8 / 0 |
| Composite context | 41 / 0 / 1 / 0 | 45 / 2 / 1 / 1 |

Full source denominators are 100 La-conditioned and 125 Ca-conditioned;
respectively 58 and 76 lack DFT pairs in this snapshot. Missing execution
results do not count as model failures or comparative accuracy improvements.

All seven DFT noncorrect cases are below. All have expected La class. `I` is
inconclusive; La/Ca denote the corresponding protocol-specific supported call.
All are seed 1. These fields uniquely identify the declared sources; exact
case IDs and unrounded DFT contrasts are in `conditioning_summary.json`.

| Protein | Conditioning | Sample | DFT | Native core | Native context | Composite core | Composite context |
|---|---|---:|---|---|---|---|---|
| A0A3F2YLY8 | Ca | 0 | I | I | La | I | La |
| A0A3F2YLY8 | Ca | 1 | I | I | La | I | Ca |
| A0A3F2YLY8 | Ca | 3 | I | I | Ca | I | I |
| A0A3F2YLY8 | La | 0 | I | I | La | I | La |
| A0ACD6B9F2 | Ca | 4 | Ca | I | La | I | Ca |
| A0ACD6B9F2 | La | 4 | Ca | I | La | I | I |
| MMOL1770 | Ca | 1 | Ca | I | Ca | I | La |

## Same completed protein summaries

- **La4:** all nine available DFT protein pools are correct in every method.
  They are A0A3F2YLY8, A0ACD6B9F2, C5ATJ3, C5AXV8, C5B120, I0JWN7,
  A8R3S4, ATQ70401.1 and BBL57595.1. The other 16/25 DFT pools are unavailable.
- **Ca5:** on nine available DFT pools, DFT, native core and composite core
  each give eight correct and one inconclusive. A0A3F2YLY8 is the inconclusive
  case; both context methods correctly call it La-supported. Native context
  gives nine correct. Context composite gives eight correct and one unavailable:
  A8R3S4 loses coverage through its existing Ca-conditioned sample 3 SCF failure.
  The other 16/25 DFT pools remain unavailable.
- **Balanced:** all eight available DFT pools are correct, as are native core,
  native context and composite core. Context composite gives seven correct and
  the same unavailable A8R3S4. The other 17/25 DFT pools remain unavailable.

Every protein pool still requires its complete declared source set. The saved
matched summary identifies all included proteins and every missing member.

## Exact same three-fold subsets

On the 37 triples presently available to DFT:

| Method | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Preserved DFT | 37 | 0 | 0 | 0 |
| Native MACE core | 37 | 0 | 0 | 0 |
| Native MACE context | 35 | 0 | 2 | 0 |
| Composite core | 33 | 0 | 4 | 0 |
| Composite context | 37 | 0 | 0 | 0 |

The full denominator remains 100: DFT has 37 correct and 63 unavailable.
Nine proteins have all four DFT triples correct; MMOL1770 contributes one
additional supported triple. Native-context inconclusives are two A8R3S4
subsets; core-composite inconclusives are two each from A0A3F2YLY8 and A8R3S4.
Original complete MACE totals remain unchanged: native core and context
composite each 94 correct/6 unavailable, native context 90/2/2/6, core composite
87/0/4/9 (correct/wrong/inconclusive/unavailable).

## Artifacts and limits

All new artifacts are under
`workspaces/accommodation_fold_DFT_20260920/partial_comparison_v2/`:
`dft_collection.json`, `fixed_comparison.json`, `triples/result.json`,
`matched_summary.json`, `conditioning_summary.json`, their command logs and
generated reports. These summaries retain all denominators, member IDs,
missing statuses, exact source
hashes, full original counts and matched counts. Existing pinned collectors
and comparators ran successfully; denominator assertions pass. Old snapshots
and running executors were not changed.

The 189 completed endpoint receipts sum to 79,186.917749 endpoint wall seconds
and 1,266,990.683984 endpoint-rank-seconds. These are partial endpoint totals,
not final allocated-node cost or a matched-hardware speed comparison.

All structures are repeats of consumed functional-class references, not new
direct affinity observations. Completion order makes this an unrepresentative
partial panel; overlapping triples are correlated. The concrete positive
finding is an A0A3F2YLY8 Ca5 coverage gain from context; on completed La4 and
three-fold summaries, preserved DFT already matches the strongest MACE arms.
Wait for the existing final collection before drawing a panel-level conclusion.
