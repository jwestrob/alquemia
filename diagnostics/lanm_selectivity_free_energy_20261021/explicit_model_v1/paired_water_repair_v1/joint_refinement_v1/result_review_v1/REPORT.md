# Bounded joint water refinement: actual result

All separate development controls pass. Stop: `all_existing_controls_passed` after 2 complete outer steps.

Objective 0.296735381 → 0.141985501; balanced exchange error -0.049294287 kcal/mol (limit 1).

| Control | Selected RMS | Selected maximum | Pass |
|---|---:|---:|---|
| La | 1.151387 | 4.838811 | True |
| Dy | 0.927327 | 3.559320 | True |
| paired | 0.276445 | 1.279527 | True |

Rigid-response units are kcal/mol/Å; torques use the unchanged 1 Å divisor. RMS/max limits remain 2/5.

## Completed steps

| Step | Selected trial | Objective | La max | Dy max | Paired max | All pass |
|---|---|---:|---:|---:|---:|---|
| 1 | bounded_LS_full | 0.215458777 | 5.624842 | 4.333684 | 1.291195 | False |
| 2 | bounded_LS_full | 0.141985501 | 4.838811 | 3.559320 | 1.279527 | True |

Passing but unselected trial slots: 1. Every trial's controls and all such rows are retained in REVIEW.json. Selection follows the original objective with current-first ties.

## Parameters and cost

| Parameter | Scout start | Selected | Move | Bound hit |
|---|---:|---:|---:|---|
| La_R | 3.686890203 | 3.710148544 | +0.023258341 | False |
| La_damping | 0.190937500 | 0.200160779 | +0.009223279 | False |
| Dy_R | 3.345904106 | 3.389554169 | +0.043650063 | False |
| Dy_damping | 0.144625000 | 0.167125000 | +0.022500000 | False |

New vectors: 24/144; native invocations: 128/1152 across 64 new cells; 158 inherited cells preserved.
Recorded wall 33.555 s; parent CPU 23.594 s; aggregate native child user/system CPU 25.564/34.005 s. Concurrent receipt CPU deltas are not summed.

All selected raw printed energy/gradient records replay exactly; sources, inputs, native receipts and parent ledgers retain their pins.

## Method sensitivity and limits

| Response | dE/dq | UnCP PBE0 residual | UnCP CC residual | CP CC residual |
|---|---:|---:|---:|---:|
| La | -4.991018 | -0.818970 | 4.774096 | 7.472010 |
| Dy | 31.664884 | -0.996566 | 3.247741 | 5.831699 |
| paired | 36.655902 | -0.177596 | -1.526355 | -1.640311 |

Positive q moves water outward; units are kcal/mol/Å. CC/CP comparisons remain separate from the objective. The CC deterioration flag is false relative to the original model, while individual CC residuals remain substantial.

This is a controls-pass stop on four consumed water fixtures, not global optimality or physical qualification. Both stored Jacobians have rank 4 at their iteration origins; the final selected point has no newly evaluated Jacobian. Protein/donor transfer, bulk hydration, electronic-state accuracy and affinity remain unqualified. Reserved validation is unopened; this report releases no further evaluation.

Result SHA256: `fab7016c4335a82dea983368abc6f23201dd83df0324f7794e5e0c978fd8e190`
