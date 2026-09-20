# Fold conditioning: solvent helps the La-source scanner, with limits

**The frozen solvent correction removes all four wrong calls among the 100 new
La-conditioned reference folds.** Context scoring changes from native MACE
93 correct/4 wrong/1 inconclusive/2 unsupported to composite
97 correct/0 wrong/1 inconclusive/2 unsupported. The declared La4 median improves
one complete protein pool: Q9Z4J7 changes from false La to correct Ca. Both methods
cover exactly the same 23/25 proteins; composite calls all 23 correctly.

This is a useful robustness result for the La-conditioned scanner workflow.
It is not a new biological validation: all 250 structures belong to the same
25 previously consumed PQQ reference proteins. The 25 canonical source replays
remain 25/25 in both core and context; the primary source-transfer set has 225
other geometries. No thresholds, labels, source selection or default changed.

## Actual primary source transfer

All rows below retain the full 225-source denominator. Counts are
**correct / wrong / inconclusive / unavailable**.

| Representation | Native MACE | MACE + frozen solvent transfer |
|---|---:|---:|
| Original core | 199 / 0 / 9 / 17 | 189 / 0 / 18 / 18 |
| Complete context | 201 / 6 / 1 / 17 | 203 / 2 / 2 / 18 |

The correction fixes five of the six native-context errors and makes the sixth
inconclusive. However, it introduces two new false-Ca calls among Ca-conditioned
folds: A0A3F2YLY8 sample1 and A0ACD6B9F2 sample4. Therefore it does not establish
conditioning-independent accuracy. On the original core it increases inconclusive
calls; solvent addition is not uniformly helpful across representations.

## Declared ensemble descriptors

These are structural descriptors using the unchanged old bands as transfer tests.
La4 is the median of all four noncanonical La-conditioned folds; Ca5 uses all five
Ca-conditioned folds; balanced is the equal mean of those two medians. Every
member must be available. They are not thermal populations or free energies.

| Context descriptor | Native MACE | Composite |
|---|---:|---:|
| La4 | 22 / 1 / 0 / 2 | **23 / 0 / 0 / 2** |
| Ca5 | 22 / 0 / 0 / 3 | 21 / 0 / 0 / 4 |
| Balanced | 21 / 0 / 0 / 4 | 20 / 0 / 0 / 5 |

The La4 gain is real within the declared transfer test: Q9Z4J7 is corrected on the
same 23-protein coverage. Ca5 and balanced already had no native errors. Their
composite coverage decreases because one additional source fails the solvent
calculation; no extra ensemble accuracy is demonstrated on their common covered
proteins. The balanced core descriptor is also 20 correct/5 unavailable.

On the same 20 complete context pools, the median absolute Ca5-versus-La4 shift
falls from 3.9134 to 2.3111 model-kcal/mol. Eleven proteins improve and nine worsen;
the largest shift grows from 7.2802 to 9.1771. Typical conditioning sensitivity
improves, but a uniform reduction is not established. A balanced average also
reduces its distance to either arm by construction; that algebra alone is not
predictive validation.

## State, failures and actual work

Preparation supports 233/250 sources. Sixteen fail the exact selected-PQQ gate;
one has overlapping expanded atoms/caps. Original failures remain in the source
record. Seven canonical cases differ from archives only by ~2e−16 Å cap arithmetic
on the CPU node; the separately recorded 1e−12 Å reconciliation changes no
coordinates or scores. Core composition/charge/cofactor/water state is invariant
on all 21 fully prepared proteins. Context membership varies in 21/25 proteins;
context differences cannot be attributed solely to geometric response.

The actual warm MACE run produced 834 fresh endpoints, including two successful
canonical execution checks, with 98 other endpoints reused. Native GFN2 attempted
1668 fresh endpoints and reused 196 exact scientific keys. Two fresh vacuum
endpoints fail native SCF after 125 cycles: P15279 La-conditioned sample3, core La;
and A8R3S4 Ca-conditioned sample3, context La. Both failures remain unavailable;
their final nonconverged energies are never used. Valid paired ALPB components
remain recorded independently. No molecular retries were made for this comparison.
The full result contains 464/500 composite case-representations; this is
232/250 in each representation, with different failed sources.

The collector verified fresh receipts and explicit reuse mappings, preserving
original task/manifest identity. Reuse required identical scientific keys and
physical input states; native inputs were checked against the current prepared
source coordinates and recorded reconciliation. All 25 canonical composite
replays exactly reproduce their archived values in both representations.

The reporting task itself launched **zero molecular calculations**. Root's
execution records under `folds_v1/` retain allocation costs and finite manifests.
Nine real-fixture parser/accounting checks passed before final collection; the
final complete-result check is recorded in the final test log. No failed or
unsupported source was removed from any denominator.

## Recommendation and artifacts

Pursue the La-conditioned ensemble as an optional robustness descriptor, while
preserving the promoted single-source mode. Do not promote a new ensemble default
or claim broad La/Ca affinity validation from this result. The source-conditioned
Ca failures and missing-pool coverage remain concrete limitations.

- [Frozen plan](FOLD_ROBUSTNESS_PLAN.md), [compact results](FOLD_COMPARISON_RESULT.json),
  [all count tables](FOLD_COMPARISON_TABLES.md), [commands](FOLD_COMPARISON_COMMANDS.md).
- Full unrounded components, all source errors and protein pools:
  `workspaces/accommodation_goal_20260920/folds_v1/comparison_v1.json`.
- Flat 500-row export: `folds_v1/comparison_export_v1/sources.csv`.
- Source, preparation and actual collections are hash-pinned in the full result.
