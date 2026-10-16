# A missed amide orientation explains one XoxF failure

**The uniform search repairs one unique XoxF source, preserves every previously
correct call, and removes the sampled Ca/La score overlap.** It does not resolve
the other eight abstentions. The evidence supports a specific local-search
limitation; it does not show that the folds as a whole are bad or that all remaining
uncertainty shares this cause.

All 257 current source/context pairs have complete scores. Their 260 cohort
memberships include three aliases. These are repeatedly examined development
structures, including the designated calibration inputs, not 257 independent
proteins or a blind validation set. The 18 superseded historical records remain
separately visible in the original 278-record ledger.

## What actually improved

| Method | Correct unique contexts | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Corrected frozen-charge three-candidate comparator | 248 | 0 | 9 | 0 |
| Add the intact fixed amide alternative | 249 | 0 | 8 | 0 |
| Also add paired bounded relaxation | 249 | 0 | 8 | 0 |

Operational and mathematical selection give the same counts. The sole repaired
source is **A0ACD6B9F2, Ca-conditioned seed1/sample4**. Its appearance in both
cohorts explains the two cohort-level improvements: expanded135 changes 129/6 to
130/5 correct/inconclusive, while currentCa125 changes 121/4 to 122/3. Those gains
must not be added as two distinct repairs. No correct call or score coverage is lost.

Here the comparator is the corrected-parameter, frozen-charge research protocol,
not a claim that these changes have been promoted to the production default.

The combined known-class raw ordering gap changes from **−1.370431 to
+3.633118 kcal/mol** already with the fixed alternative. Subsequent relaxation
does not increase that gap. This is a descriptive comparison of consumed data;
no threshold was chosen from these extrema. All original canonical bands and
the **8.450323470286094 kcal/mol** canonical gap remain unchanged, as do all
25 canonical scores and selected geometries.

The existing 100 correlated three-source groups retain **89 correct, zero wrong,
five inconclusive and six unavailable** for every tested pool. Do not replace
this result with an older, differently parameterized all-correct ensemble claim.

## Why the XoxF source was different

The source presents the terminal Asn amide nitrogen toward the metal and its
oxygen away. The original bounded torsion search could not reach the intact
alternative orientation. Offering that alternative to **both** metals adds a
competitive local arrangement without changing the atoms, proton inventory,
cofactor state, solvent inventory, source context or thresholds.

This is not proof of atom misassignment or a defective whole-protein fold.
Nor does it identify Asn orientation alone as the entire energetic cause:
the subsequent bounded search moves other donors too. The fixed alternative
retains a short represented hydrogen contact. Relaxation lengthens the named
contact, retains the oxygen-facing orientation, and improves source consistency.

Only A0AC Ca4 and La4 change scores anywhere in the full transfer. From the fixed
alternative to the relaxed pool, the A0AC Ca-conditioned range contracts from
12.350985 to 5.677901 kcal/mol; its La-conditioned range contracts from 19.786377
to 2.222541. These are matched-source ranges, with no lost coverage. The extra
relaxation therefore provides a real structural-consistency benefit, despite
adding no further correct class call.

Independent four-cell native DFT supports the large La stabilization and net
direction on the selected relaxed paths: Ca work +1.906580, La −17.400253,
paired change +19.306832 kcal/mol, versus MACE −0.286036, −18.274380 and
+17.988344. DFT contradicts the small Ca stabilization. This is conditional
local-energy support, not exact agreement, validation of the unrelaxed fixed
flip, solvent accuracy, a DFT pool minimum, electronic stability or an affinity.
See [the reference result](../reference_parallel_v2/FINAL_REPORT.md).

## What remains uncertain

The eight remaining current abstentions are:

- PedE Q88JH5: La-conditioned samples3 and4.
- ExaA Q9Z4J7: La-conditioned samples2,3,4 and two-Ca sample3.
- MMOL1770: Ca-conditioned sample0.
- C5AXV8: two-Ca sample4.

Their scores now remain correctly ordered against the other known class in this
development collection, but fall inside the unchanged conservative canonical
interval. That explains their immediate classification status; it does not identify
a unique physical cause of their score shifts. Shrinking that interval after
seeing these labels would not establish new predictive accuracy.

The completed case-specific work remains relevant. The PedE experiments expose
a solvent-response/calibration tradeoff; they do not prove that all electronic
response is spurious. MMOL and the tested ExaA motions retain their response under
independent DFT, weakening a MACE-only explanation. The tested packing, reciprocal
metal, midpoint, small PQQ-pose and cavity remedies did not supply a useful rescue.
Those results do not rule out every environmental or chemical-state effect.
See [the causal synthesis and its limits](../../review/REPORT.md).

## Practical decision

Retain the chemically meaningful alternative-start mechanism as a supported
research improvement. Do not promote the unrelaxed fixed flip merely because its
class counts match the more expensive search: its contact geometry and its distinct
reference-validation status still matter. Uniform relaxation has not demonstrated
additional classification benefit beyond its A0AC source-consistency gain.

A separate **retrospective** policy replay is complete. The existing composite
fixed-flip selection screens the alternative without using class labels. If either
metal selects it, both bounded searches and shared-pool cross-scoring are enabled;
otherwise the original three-candidate result is retained. Missing required
results remain unavailable. Only A0AC Ca4 and La4 trigger among257 contexts.
The routed result exactly reproduces full common6 scores, both metal selections,
decisions and availability under both reporting policies. Thus the geometry and
source-consistency benefit can be retained by this rule on the consumed archive.

This supports a live opt-in screen-and-relax candidate, not promotion of the
short-contact fixed geometry. Four eligible and510 avoidable search slots are
retrospective counts, not measured savings. The screen can miss an initially
unfavorable alternative that becomes favorable after relaxation. No prospective
accuracy, live runtime improvement or production promotion is established.
See [the rule, result and limitations](gated_replay/REPORT.md).

## Executed work and cost

| Stage | Job | Scheduler seconds | CPUs | GPU |
|---|---:|---:|---:|---|
| New-context MACE scout | 5263 | 193 | 12 | 1 A5000 |
| Remaining MACE search | 5297 | 9404 | 12 | 1 A5000 |
| Solvent preflight | 5610 | 1610 | 7 | none |
| Full new solvent components | 5661 | 24381 | 7 | none |
| Saved report | 5666 | 203, cancelled after application completion | 1 | none |
| Editable figure | 5748 | 13 | 1 | none |
| Saved deduplicated summary | 6035 | 1 | 1 | none |
| Saved gated-policy replay | 6045 | 2 | 1 | none |

The solvent stage admits all **1338 new unique endpoints**, with **5352 SCC calls
and 2676 native component calls**, reusing 1746 compatible endpoint results. The
remaining MACE stage records 7071 model requests and no failed requests; its actual
forward count is unavailable. Exact internal timings and source receipts are in
[the compact result](delivery/RESULT.json). Scheduler wall time differs slightly
from application timing. These are batch development measurements, not a matched
production latency comparison or measured CPU/GPU utilization. Preparation,
preflight, original failures and reused campaign costs remain in their original
receipts; this table is not a claim of total campaign cost.

Reporter5666 finished its application in 4.19 seconds, then its allocation
lingered. Only that completed owned allocation was cancelled; all report artifacts
remain valid. Its scheduler cancellation is preserved separately. Figure5748 and
summary6035 subsequently completed normally. No molecular calculation was repeated.
Replay6045 also completed normally, with zero molecular calls. Its2 allocated
CPU-seconds measure the saved-data replay only, not the screened molecular workflow.

## Reproduction and artifacts

- [All individual contexts, PNG](delivery/COMMON3_COMMON6_operational.png),
  [editable SVG](delivery/COMMON3_COMMON6_operational.svg),
  [PDF](delivery/COMMON3_COMMON6_operational.pdf),
  [figure data](delivery/FIGURE_DATA.tsv). The actual PNG was visually inspected.
- [Compact results and exact source pins](delivery/RESULT.json).
- [Saved-data decision summary command](decision/READY.md); this performs no new
  molecular calculation. Use a new output directory for a deliberate replay.
- [Original staged preparation/search/solvent commands](READY.md). These are
  already executed; their availability is not a request to repeat them.

Authoritative solvent collection SHA256:
`09ee0ae6f26b8e7f479d695c6c4b51da49bbf5a656793d1dbe54e93300c994a1`.
Deduplicated decision result SHA256:
`f70378b775db2844a18f1cbfb5a974b3b3eb6b930aab2b0ea65a0e6831ba61eb`.
Full outputs remain under `workspaces/pqq_failure_causal_synthesis_20261015/` in
`asn_transfer_solvent_v1`, `asn_transfer_report_v1`, `asn_transfer_figure_v1` and
`asn_transfer_decision_v1`. Production defaults and immutable earlier results are
unchanged. No new calibration or full-cohort production rescore occurred.
