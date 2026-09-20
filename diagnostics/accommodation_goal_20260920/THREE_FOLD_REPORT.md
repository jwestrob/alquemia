# Three-fold addendum: context composite is stable, and native core matches it

**All 94 available three-fold context-composite descriptors classify correctly;
the simpler native-core descriptor also achieves 94/94 on the same coverage.**
The full denominator is 100 declared triples, with six unavailable triples in
each of those methods. This is robustness to subset choice on consumed reference
structures, not proof of superiority over native core or validation of the actual
three-fold PLM inputs.

This addendum was specified after the original frozen fold comparison. It
computes every three-of-four subset of each protein's four noncanonical
La-conditioned samples: four triples per protein, 25 proteins, both representations
and both methods. No favorable triple was selected. Every triple requires all
three scores, and every method uses its original exact bands without refitting.
There were **zero new molecular calls** and no production/default changes.

| Representation / method | Correct | Wrong | Inconclusive | Unavailable | Proteins with all four triples correct |
|---|---:|---:|---:|---:|---:|
| Native core | 94 | 0 | 0 | 6 | 23/25 |
| Core + solvent | 87 | 0 | 4 | 9 | 20/25 |
| Native context | 90 | 2 | 2 | 6 | 21/25 |
| Context + solvent | **94** | **0** | **0** | **6** | **23/25** |

The context solvent correction fixes Q9Z4J7's two wrong triples and A8R3S4's two
inconclusive triples. Each of those proteins has all four composite triples
correct. It improves native-context subset robustness; it does not outperform the
simpler native-core comparator on this panel. Adding solvent to the original core
is less useful here and remains separately reported.

## Missing coverage and worst outcomes

MMOL1770 and Q88JH5 each have one unsupported La source. For each protein, three
of four triples therefore remain unavailable; the one fully supported triple is
correct. Neither protein is counted among the 23 with every triple correct, and
no failed source was removed to produce a selected success. Core + solvent also
inherits the P15279 SCF failure, adding three unavailable triples. All individual
triples, omitted members, scores and missing members are retained in the export.

Per-protein summaries retain all flags/counts, a worst-scored outcome and a worst
outcome including missing coverage. The latter uses the predeclared order wrong,
unavailable, inconclusive, correct. Correlated subsets are not 100 independent
biological observations. This analysis reuses already-inspected energies and is
explicitly developmental.

## Use and limits

A strict three-fold median is a reasonable opt-in candidate for the existing
La-conditioned scanner: on these references, choosing any supported triple does
not damage context-composite calls. Native-core parity means this result alone
does not justify the extra context/solvent cost over every simpler option.
The actual PLM three-sample policy still needs its own input compatibility and
coverage checks; PLM predictions remain unlabeled. No ensemble promotion or
rescore is performed here.

- [Recorded rule](THREE_FOLD_PLAN.md), [compact counts and per-protein worst outcomes](THREE_FOLD_RESULT.json),
  [every triple/method/score](THREE_FOLD_ALL_TRIPLES.csv), [four real-fixture tests](THREE_FOLD_TESTS.txt).
- Full immutable result: `workspaces/accommodation_goal_20260920/folds_v1/triples_v1/result.json`.
- Original La4/Ca5/balanced results remain unchanged in
  [the frozen report](FOLD_COMPARISON_REPORT.md).

Reproduce into a new output file from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 scripts/accommodation_folds_triples.py \
 --comparison workspaces/accommodation_goal_20260920/folds_v1/comparison_v1.json \
 --plan diagnostics/accommodation_goal_20260920/THREE_FOLD_PLAN.md \
 --output workspaces/accommodation_goal_20260920/folds_v1/triples_v1/replay.json
```
