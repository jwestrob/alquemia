# Alternative amide search retains canonical discrimination

The complete solvent calculation preserves all 25 canonical reference calls.
Adding the fixed terminal-amide alternative and its two bounded relaxed proposals
changes no canonical selected geometry or paired contrast. The operational and
mathematical calibration gaps remain exactly 8.450323470286094 kcal/mol; the
designated-reference bands equal the original bands. No transfer outcome was
used to move a threshold.

| Population | Original common3 | Expanded common6 |
|---|---|---|
| Designated canonical references | 25 correct, 0 wrong, 0 inconclusive | Unchanged |
| A0ACD6B9F2, ten consumed source folds | 9 correct, 1 inconclusive | 10 correct |

The A0AC result is a structural robustness gain within one protein, not ten
independent biological validations. Its Ca-conditioned score range contracts
from 23.326 to 5.678 kcal/mol, and its La-conditioned range from 20.459 to 2.223.
Both selection policies agree. The separate origin-only comparator remains
8 correct, 1 wrong and 1 inconclusive; it is not silently replaced by the search.

## What was executed

Canonical solvent job5158 completed 144 missing endpoint geometries: 576 SCC
backend calls and 288 native frozen-charge component evaluations. It reused
150 original cells and six exact A0AC additional cells. All 300 canonical matrix
cells pass their declared component/native admission. No new MACE, DFT or
optimization ran in this solvent job. Scheduler cost was 2644 seconds on seven
CPUs (18508 allocated CPU-seconds), zero GPUs. The internal measured duration
was 2643.268097995897 seconds. Original failed preflight5135 and successful
recovery5157 remain preserved separately.

Saved-summary job5202 independently checks the actual failed5135 missing-result
fixture, replays the completed5030 A0AC scores/selections/ranges, and verifies
the canonical/A0AC shared source's exact state, mappings and native receipts.
Both real-artifact checks pass. It performs no molecular evaluations or refit.

Scientific acceptance uses the already qualified component-v2 comparison and
native four-component replay. Historical absolute-v1 parity failures remain
failures; they were not reclassified or used as the current admission rule.
Boundary flags on unselected proposals remain recorded. These are finite
constrained search results, not unconstrained minima or equilibrium populations.

## Decision and next work

The canonical gate permits the declared uniform transfer test. Apply the same
anchor-Asn alternative to the complete current expanded135 and corrected Ca125
cohorts, retaining 260 memberships across 257 distinct source/context pairs.
Reuse the 34 already completed unique contexts; preserve all exclusions and
failed additional cells. This tests whether the local improvement transfers
without losing correct Ca or La calls. Production is unchanged.

Exact staged commands and source pins: [transfer readiness](../transfer/READY.md).

## Authoritative artifacts

- `workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_solvent_v2/execution/FINAL_COLLECTION.json`
- The neighboring `EXECUTION.json` records all new calls and resource receipts.
- `workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_solvent_report_v1/execution/RESULT.json`, SHA256 `45f5f8751ac1b123cc6b042162987492fa699ca4ccb7db1aa4c0abefd01e9a10`.
- Neighboring `CANONICAL.tsv`, `A0AC.tsv`, and `A0AC_RANGES.tsv` contain editable results.

All sources are consumed development evidence. Canonical accuracy alone does
not establish broader biological transfer.
