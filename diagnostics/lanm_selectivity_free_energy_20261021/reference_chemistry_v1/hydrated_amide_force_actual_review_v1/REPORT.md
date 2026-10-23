# Hydrated acetamide: small steps support the analytic gradient

**Both metals and their difference pass the frozen numerical checks.** The original La anomaly at ±0.15 Å is consistent with finite-step curvature along this ligand translation. At ±0.010/±0.005 Å, the largest total-energy finite-difference versus analytic discrepancy is **0.00369 kcal/mol/Å**, and the largest two-step change is **0.00483 kcal/mol/Å**. Both are below the declared 0.1 resolution.

The old La coarse-screen failure remains preserved: analytic 14.852538 versus secant 13.139391, a difference of 1.713148 exceeding its original 1.485254 screen. Its archived status and thresholds have not changed.

## Actual derivative comparison

Derivatives are in kcal/mol/Å. Positive motion translates the intact nine-atom ligand outward while the metal and seven waters stay fixed. Analytic central derivatives are La **14.852538**, Dy **31.729427**, and Dy−La **16.876889**.

| Quantity | Step Å | Electronic | D4 | Total | Total FD − analytic |
|---|---:|---:|---:|---:|---:|
| La | 0.010 | 13.017099 | 1.834304 | 14.851396 | −0.001143 |
| La | 0.005 | 13.021920 | 1.834273 | 14.856225 | +0.003686 |
| Dy | 0.010 | 30.059098 | 1.669458 | 31.728559 | −0.000868 |
| Dy | 0.005 | 30.062311 | 1.669489 | 31.731788 | +0.002361 |
| Dy−La | 0.010 | 17.041998 | −0.164847 | 16.877163 | +0.000274 |
| Dy−La | 0.005 | 17.040391 | −0.164784 | 16.875564 | −0.001325 |

Electronic, D4 and total two-step checks all pass separately. The paired total two-step difference is −0.001600. The 0.005 Å results therefore qualify as the prespecified summaries. Small component-closure differences reflect printed native energy precision and remain within the unchanged native admission criterion.

## Evidence and limits

All eight native cells are admitted. The independent review passed **122 checks**, replaying the original collector exactly, checking actual input/geometry/job joins and recomputing derivatives from printed energy tokens with decimal arithmetic. No COSX/final-integration stage signatures changed between the saved central EnGrad calculations and displaced single points; the expected presence/absence of the gradient stage remains explicit.

This qualifies the sampled central ligand-translation derivative at the declared resolution. It does not establish all-coordinate force accuracy, physical electronic-state accuracy, basis convergence, or eligibility of other donor states. It neither validates the force-field parameters nor changes the separate failed transfer result.

## Actual cost and pins

Jobs 17579–17586 used eight single points, four CPUs each, **14,016 allocated CPU-seconds (3.89 CPU-hours)**. Their measured native elapsed-time × allocation sum is 13,987.37 CPU-seconds. Collector17587 used another **2 allocated CPU-seconds**. There were no GPU calls or repeated centers. These allocation products are not measurements of CPU utilization. This saved review took 4.73 seconds on one local thread and made no native calls or submissions.

| Artifact | SHA256 |
|---|---|
| Actual review RESULT.json | `a6128e37bdd14d6d9c6ae4d93c5cde394bf6e971ed50b00b75f3d70c4194c7e6` |
| Review source | `b55a748fe05ccc3a27fd7cb823d7661170f506f960ea03a46fdb1b00cbd37d14` |
| Original collector RESULT.json | `ac407b50fc1b278b44906e69b35fc25a32734599430f90228f8fb7d875e3c241` |
| Frozen native adapter | `7de9c17aa9198ba07c66c3f934390fc2a67d94a77860eca8ba4c0378b25ad677` |
| Frozen eight-cell PANEL.json | `38b42c1020c47b12512076d59d66e299d4c8eb7f14754f8cfe5f71bec502e4cc` |

Review result: `workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_amide_force_actual_review_v1/review_v1/RESULT.json`. Its raw-input, output, native residual, energy-stage and scheduler-receipt pins preserve the complete evidence.
