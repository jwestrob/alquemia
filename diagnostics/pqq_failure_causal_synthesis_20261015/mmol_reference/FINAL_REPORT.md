# MMOL1770: independent reference reproduces Ca-directed accommodation

All **6/6 declared vacuum endpoints** passed the frozen energy-diagnostic collector. Native ORCA6.1.1 r2SCAN-3c reproduces the Ca-directed accommodation of MMOL1770 Ca-conditioned sample0. A MACE-only artifact is therefore unsupported by this comparison; quantitative agreement is incomplete.

The exact 193-atom/13-cap dry oxidized-PQQ3− context, Ca−3/La−2 singlets and saved origin/adaptive_Ca/adaptive_La geometries were retained under the [frozen plan](../MMOL_REFERENCE_PLAN.md). No DFT optimization or solvent calculation occurred; source1 remains inventory-only.

## Actual response

Work is relative to each metal's own origin; contrasts are Ca minus La, in kcal/mol.

| Comparison | Native DFT | MACE | DFT − MACE |
|---|---:|---:|---:|
| Same adaptive_Ca geometry | −10.822777732 | −9.136199735 | −1.686577996 |
| Same adaptive_La geometry | +1.216792487 | +1.989215465 | −0.772422978 |
| Complete common3 vacuum pool | −5.221984827 | −4.004941242 | −1.217043585 |

Both methods select **Ca@adaptive_Ca and La@adaptive_La**. Their selected Ca/La works are −14.424446311/−9.202461484 (DFT) and −14.319518665/−10.314577422 (MACE). These mixed pool minima are distinct from either same-geometry path. The historical composite used the same selections, but its −4.186093026 response includes solvent. The Ca path and pool discrepancies exceed the declared 1 kcal/mol interpretive resolution; this resolution is neither an error bar nor a fitted acceptance threshold.

## Electronic-state limits

Last endpoint **4639, La@adaptive_Ca**, terminated normally after79 cycles. Its final TRAH/NR residual 5.899231e−6 passes the printed 1e−5 criterion. Energy change −1.2178e−9 passes1e−8, but MAX/RMS density changes .089631/.00059696 exceed their printed tolerances. Its final printed frontier ordering is occupied419 −.026826Eh, empty420 +.073174Eh; this does not establish electronic stability. Other endpoints retain energy/gradient-based endings, large density residuals and, in several cases, inverted final frontier ordering. The inherited collector accepts energy diagnostics; same-state continuity, stability, force accuracy and a unique ground state remain unqualified.

## Actual allocation costs

`sacct --allocations`, 2026-10-15 03:55:36UTC; all COMPLETED/0:0. CPU-seconds below are allocated, not measured process CPU.

| Job / endpoint | Wall seconds | CPUs | Allocated CPU-seconds |
|---|---:|---:|---:|
|4603 Ca origin|3978|12|47736|
|4604 La origin|3789|12|45468|
|4637 Ca adaptive_Ca|4093|12|49116|
|4638 Ca adaptive_La|4585|12|55020|
|4639 La adaptive_Ca|5497|12|65964|
|4640 La adaptive_La|3683|12|44196|
|4641 comparison|0|1|0|

## Decision and provenance

Prioritize the shared source/context/state assumptions before attributing the inconclusive composite result to MACE alone. This finite comparison establishes neither affinity nor classifier improvement and does not validate DFT as an oracle. No automatic follow-up is launched.

Authoritative [RESULT.json](../../../workspaces/pqq_failure_causal_synthesis_20261015/mmol_comparison_v1/execution/RESULT.json), SHA256 `cf0216391de0680b50b746c9799edd5365c23cc6ee1b9855e39844173903ab61`, retains all six input/output/receipt pins and unrounded values. Final4639 output SHA256 `241f99b87f9f87f4d40dfac53f7cbbdae33cecbf4f7887efc8100b415087bf41`. This report reads saved outputs and scheduler metadata only.
