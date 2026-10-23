# The model's electrical-response errors partly cancel

The unchanged water model underestimates the separate metal–water dipole response and the collective change when those waters form a shell. These errors partly cancel in the full-shell response. The collective La/Dy ordering also reverses relative to this quantum reference.

C is complex minus ion minus intact water shell; P is the sum of individual metal–water interaction differences; B = C − P. All are neutral dipole differences in the same parent coordinate frame. The table gives the z component in native Tinker Debye units; transverse components are negligible in these fixtures.

| Water shell | Quantity | Quantum | Model | Model − quantum |
| --- | --- | ---: | ---: | ---: |
| W9 La | Full interaction C | +2.244444 | +2.107200 | −0.137244 |
| W9 La | Pair sum P | +3.929194 | +2.916725 | −1.012469 |
| W9 La | Collective difference B | −1.684749 | −0.809525 | +0.875224 |
| W9 Dy | Full interaction C | +2.139823 | +2.017552 | −0.122271 |
| W9 Dy | Pair sum P | +3.959697 | +2.729250 | −1.230447 |
| W9 Dy | Collective difference B | −1.819874 | −0.711698 | +1.108176 |

The full interaction is about 6% too small, while the collective magnitude is about 52% too small for La and 61% for Dy. For this W9 geometry, Dy−La in B is −0.135125 Debye in the quantum result and +0.097827 in the model. These are electrical moments, not energy differences or binding preferences.

All four parents are retained. W8 model vectors cancel by symmetry; quantum C/P/B z values are La +0.000965549/+0.000124167/+0.000841382 and Dy +0.000832664/+0.000045980/+0.000786684 Debye. Some exceed print bounds, but without a quantum numerical-error bound they do not establish physical symmetry breaking or useful model accuracy.

## Admission and precision

All 78 model cells are available: four successful native outputs from job18313 were recovered after filename/type-label reader failures, and job18315 computed only the remaining74. Original source, failed application result and native outputs remain preserved. Both bare ions passed before the other72 new calls. Independent actual review passes4,747 checks; all38 saved energies replay exactly. Maximum induction residual is4.3e−9 Debye under the unchanged1e−8 limit. Permanent contributions cancel atom-by-atom with identical local frames; rounded total moments provide a separate consistency check.

The quantum reconstruction and independent full-precision review admit all80 records. Shared parent coordinates compare exactly. Conversion uses the pinned native constants,2.541747577148327 Debye per atomic dipole unit. Combined print bounds per W9 component are3.563e−5(C),1.047e−4(P),1.317e−4(B) Debye, far below the W9 discrepancies. These bounds cover serialization only, not SCF or reference-method error.

Native scheduler costs:18313,40 allocated CPU-seconds;18315,92; total132. The comparison itself used32.350 seconds wall/31.308 CPU-seconds on one local thread and made zero molecular calls. Independent native review used7.068 seconds wall/5.926 CPU-seconds. Reused quantum calculation costs are not charged again.

## Meaning for the next mechanics change

This is direct evidence that the model misses electrical redistribution within the chosen PBE0/ECP reference. Matching the full-shell moment alone would conceal the pair and collective errors. A repair must address both contributions while retaining energy and force tests; another fit to one coordination contrast would not establish that.

It does not uniquely identify charge transfer, polarization, charge penetration or Pauli effects, nor establish that this dipole error causes the entire energy discrepancy. The population analyses disagree on the small paired charge-transfer sign. Frozen-f, finite-basis and SCF limitations remain. No parameter, target, validation set or production setting changed. A separate bounded correlated-method comparison is being prepared to test reference uncertainty before any stronger physical attribution.

Exact inputs, all Cartesian vectors, bounds and costs are in the matching workspace `actual_v1/RESULT.json`. Source review is `independent_source_review_v1/REVIEW.json`; independent actual model review is under `explicit_model_v1/saved_redistribution_counterpart_v1/independent_actual_review_v1/actual_v1`.
