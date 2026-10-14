# Midpoint gains largely cancel in the target's metal contrast

**The interior probe finds lower composite energies, but barely changes the problematic target score. Close this finite fraction test without finer sampling or joint-optimizer development.** The result establishes that the frozen metal coordinate excludes some favorable candidates; it does not establish a useful discrimination repair.

## Complete common-pool result

All eight new endpoints are available and all20 prior endpoints are reused. The historical three-candidate and previous five-candidate pools reproduce their old selections. Mathematical and operational policies give the same selections in both sources. Every work below is relative to the same metal's unchanged original origin, in kcal/mol.

| Host | Metal | Historical3/previous5 choice | New7 choice | Old work | New work | New−old |
|---|---|---|---|---:|---:|---:|
| Sample0 control | Ca | adaptive_Ca | adaptive_Ca | −7.624675 | −7.624675 | 0 |
| Sample0 control | La | adaptive_La | midpoint_adaptive_La | −12.934776 | −14.552547 | −1.617772 |
| Sample3 target | Ca | adaptive_Ca | midpoint_adaptive_Ca | −9.698819 | −11.141071 | −1.442252 |
| Sample3 target | La | adaptive_La | midpoint_adaptive_Ca | −26.390560 | −27.655842 | −1.265282 |

The target now selects the same physical midpoint/donor geometry for both metals. Their energy gains nearly cancel: the Ca−La scalar changes only **−0.176971**. The control scalar changes **+1.617772**, because only La benefits. The stored target-minus-control scalar change is **−1.794742**. Thus reduced source separation arises mainly from the control moving in the La-directed direction, rather than removal of the target's large differential accommodation. No classification is assigned; no bands are applied.

## What drives the selected changes

For target Ca, the selected change comprises MACE **−2.642239** and frozen solvent **+1.199986**. For target La, these are **−3.423274** and **+2.157992**. The target differential change is consequently MACE **+0.781035**, solvent **−0.958006**, total **−0.176971**. The solvent term offsets much of each metal's vacuum gain.

The control La choice is especially informative: its MACE energy becomes **+0.386260** less favorable while its frozen-solvent term improves **−2.004032**, yielding the net **−1.617772**. A useful composite candidate need not lower the vacuum MACE energy relative to the old selected geometry. The saved vacuum-force reversal motivated this point; it does not supply a solvent-composite gradient or establish that a vacuum-only optimizer will find the best composite candidates.

These findings retain the earlier causal metal-position sensitivity while limiting its practical interpretation. Both metals can accommodate the target's midpoint similarly, and the tested source difference remains. Lower energy and a smaller two-source spread do not demonstrate improved biological specificity, robustness across the other folds or a correct equilibrium ion position. The midpoint is a modeled, source-conditioned point with unchanged dry oxidized-PQQ chemistry; there is no independent DFT test of these new endpoints.

## Decision and qualification

End this declared one-fraction branch. Do not infer that all metal relaxation is useless, but do not escalate to finer fractions, wider bounds, a joint optimizer or recalibration on the strength of these two consumed examples. Any subsequent work needs a separate physical question and a comparison that measures metal discrimination rather than only lower energies. The original full-transfer zero-gain result remains intact.

Every new endpoint passes its recorded component-v2 acceptance and the existing native frozen-charge qualification route. The collector records156 saved-algebra checks, maximum residual **2.7436×10⁻¹¹ kcal/mol** against the unchanged10⁻⁸ tolerance. There are no missing new cells. These are numerical/implementation checks, not validation of the cofactor-state assumption or an affinity prediction.

Job4137 records **277.8259 s on seven allocated CPUs**, **1944.7815 allocated CPU-s**, eight MACE forwards,32 started/terminal SCC calls and16 native component calls. No DFT or optimization occurred. These are the current stage's costs, not end-to-end new-site throughput; underlying measured resource receipts are pinned separately.

Authoritative result: `workspaces/exaa_metal_midpoint_20261014/scoring_v1/execution/FINAL_COLLECTION.json`, SHA256 `e7490cd21dd75ce81cde99660f48ecc541a2bffc6b39d5ba7aada2388116f6e7`; cost record: adjacent `EXECUTION.json`. Config SHA256: `f0edcf081a019a5b5bffb5278050fb1ce6052361837064527a339d213b2c7147`.

All numerical values above were read from actual saved outputs. This independent review performs no scientific recomputation, new evaluation, test, submission, refit or shared-code change.
