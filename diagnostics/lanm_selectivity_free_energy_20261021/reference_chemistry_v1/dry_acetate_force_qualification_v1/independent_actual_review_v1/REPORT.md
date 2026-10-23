# Independent actual dry mono-acetate derivative review

**Passed: all 16 displaced outputs and 18 prescribed numerical checks.** Both original central gradients replay from `.engrad`; all four scout outputs are reused exactly. No center rerun or omitted cell.

| Direction | Comparison | Analytic | FD 0.010 | FD 0.005 |
|---|---|---:|---:|---:|
| radial | La | 99.58017697 | 99.55708004 | 99.55706517 |
| radial | Dy | 105.90122134 | 105.89751697 | 105.91895131 |
| radial | paired | 6.32104437 | 6.34043693 | 6.36188614 |
| rotation | La | 28.34118343 | 28.32497683 | 28.35794034 |
| rotation | Dy | 26.60499886 | 26.57697982 | 26.60124332 |
| rotation | paired | -1.73618458 | -1.74799701 | -1.75669702 |

Units are kcal/mol/Angstrom, including the fixed 1 Angstrom rotational coordinate. Maximum FD/analytic discrepancy is **0.04084177**; maximum two-step difference is **0.03296351**, both below the unchanged 0.1 resolution. Independent raw final/electronic/D4 tokens, gradient signs, torque projections, geometry transformations, input/state lineage and all component derivatives agree with the saved result.

Every displaced calculation uses Energy Check stopping. Printed density/DIIS failures remain; COSX stage signatures and DFT/XC settings agree with centers, but coordinate-dependent grid counts and central EnGrad versus displaced SP remain explicit. This qualifies only the sampled directions and resolution, not physical reference accuracy, other donor forces or a unique electronic state.

Sixteen native calls: 2464.288528 summed wall seconds and **9857.154114 native allocated CPU-s**. Saved scheduler totals: **9916 CPU-s**, plus **5 CPU-s** for the two collectors. Historical center costs are not charged again. This review made zero native calls and used one local thread.
