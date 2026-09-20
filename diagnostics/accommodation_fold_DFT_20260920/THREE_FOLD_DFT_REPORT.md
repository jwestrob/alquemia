# Preserved DFT on the same three-fold challenge

2026-09-20. The comparison is ready; the DFT campaign is still running.
This tests whether the three-fold MACE result improves on the preserved DFT
discriminator using **exactly the same sources and subsets**.

The [addition was declared before evaluation](THREE_FOLD_DFT_PLAN.md). It uses
all four three-of-four noncanonical La-conditioned subsets for each of the
25 consumed reference proteins: 100 overlapping triples. Each median requires
all three source results. The original DFT bands and four existing MACE arms
are unchanged; there is no fitting, new molecular calculation or promotion.

## Actual partial result

The immutable snapshot contains **49 complete and 367 pending DFT endpoints**,
with no recorded failures. Eighteen completed La-source endpoint receipts,
outputs and input mappings were independently verified for this comparison.

| Method | Correct | Wrong | Inconclusive | Unavailable | Proteins with all four correct |
|---|---:|---:|---:|---:|---:|
| Preserved DFT | 8 | 0 | 0 | 92 | 2/25 |
| Native MACE core | 94 | 0 | 0 | 6 | 23/25 |
| Native MACE context | 90 | 2 | 2 | 6 | 21/25 |
| Composite core | 87 | 0 | 4 | 9 | 20/25 |
| Composite context | 94 | 0 | 0 | 6 | 23/25 |

Only A0A3F2YLY8 and A0ACD6B9F2 currently have complete DFT triples in this
snapshot. Availability depends on execution order: **these partial counts do
not establish which method is better**. All MACE counts reproduce the existing
completed addendum exactly. Structural repeats and overlapping triples are
not independent biological observations or direct affinity measurements.

Artifacts under `workspaces/accommodation_fold_DFT_20260920/triples_v1/`:

- `dft_partial_collection.json`: immutable real collection snapshot.
- `comparison_with_units/result.json`: all 100 subsets, 25 protein summaries,
  exact source memberships, medians, decisions and missing-member statuses.
- `comparison_with_units/all_triples.csv`: long-form results, explicitly
  distinguishing DFT kcal/mol from MACE model kcal/mol.
- `implementation/`: pinned comparator and its five imported dependencies.
- `comparison/`: earlier identical numerical comparison, preserved before
  clarifying the CSV unit label; its source is retained separately.

## Validation and final command

Five real-artifact tests pass in 1.433 s, with zero skips
([log](THREE_FOLD_TESTS_WITH_UNITS.txt)). They check exact membership, unchanged
MACE values, real DFT energy algebra, complete-member medians, actual prelaunch
unavailability, and rejection of explicitly corrupted real bands/membership.
No scientific executable was invoked by this comparator.

Once the existing DFT collectors create `run_v1/final_collection.json`, run
from the repository root:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/accommodation_fold_DFT_20260920/triples_v1/implementation/accommodation_fold_dft_triples.py \
  --dft-collection workspaces/accommodation_fold_DFT_20260920/run_v1/final_collection.json \
  --mace-triples workspaces/accommodation_goal_20260920/folds_v1/triples_v1/result.json \
  --plan diagnostics/accommodation_fold_DFT_20260920/THREE_FOLD_DFT_PLAN.md \
  --output workspaces/accommodation_fold_DFT_20260920/triples_final_v1
```

The command creates a new immutable directory with JSON, Markdown and CSV.
Failed endpoints remain unavailable and retain their status. Existing
La4/Ca5/balanced comparisons, running snapshots and production defaults are
unchanged. The four DFT jobs continue; no extra polling or launches are needed
for this addendum.
