# Model 18473: actual conditional susceptibility admitted

All **19/19 original-baseline calls** completed, with no native or application failure. Saved-only admission rechecked exact inputs, binary/receipt lineage, all raw energies and site moments, initial controls and full numerical summary. No source, field, tolerance or output was changed and no new native call was made.

Application RESULT SHA256: `109ca5487155fe8c0e725ed59af0432adc216cc03b6bcba47bb21fd4a4c97392`.
Initial actual REVIEW SHA256: `924e1076260cc35721df6006272ca58a5d64428e4dd0f6f70d57b9d3f5deaef2` in the matching workspace `initial_actual_review_v1/REVIEW.json`. This author-side admission precedes Planck's independent actual review.

## Numerical controls

- Explicit zero-field Dy trimer: energy and every saved site dipole reproduce the baseline **printed-exactly**.
- Bare La at ±0.001 au: expected energy −0.003514119187 kcal/mol; each actual energy error is −8.13e−10 kcal/mol. Expected induced z dipoles are ±0.028468108603 Debye; maximum error is 1.086e−7 Debye. Sign, nonzero response, field conversion and the half-energy factor pass.
- All **24 state checks and 12 increment/paired checks** pass at the unchanged prospective thresholds. All 38 induction solutions admit; maximum residual is 5.9e−9 Debye.
- Largest state energy-curvature versus moment-slope discrepancy: 0.001472 a0³. Largest corresponding increment discrepancy: 0.002133 a0³. Every discrepancy lies within its combined propagated printing bound.
- Maximum state half-step energy change is 6.37e−5 a0³; maximum increment/paired energy-step change is 1.434e−4 a0³. Maximum state odd energy slope versus reconstructed zero-field moment discrepancy is 9.15e−7 ea0, below 1e−4.

These observations are consistent with the baseline's fixed-geometry linear response. They establish numerical readback/field consistency, not physical reference agreement or a general stability proof.

## Actual directional response

Values are chi_zz in a0³. Both raw field scales and both estimators are retained; no averaging or extrapolation replaces them.

| State or increment | Energy, h=.001 | Energy, h/2=.0005 | Moment, h | Moment, h/2 |
|---|---:|---:|---:|---:|
| La dimer | 18.819716 | 18.819732 | 18.819729 | 18.819729 |
| La far trimer | 26.277563 | 26.277532 | 26.277393 | 26.278377 |
| Dy dimer | 12.926084 | 12.926020 | 12.925949 | 12.925359 |
| Dy far trimer | 20.785375 | 20.785407 | 20.786092 | 20.786879 |
| La: trimer−dimer | 7.457848 | 7.457800 | 7.457664 | 7.458648 |
| Dy: trimer−dimer | 7.859292 | 7.859387 | 7.860143 | 7.861520 |
| Dy−La increment | +0.401444 | +0.401588 | +0.402479 | +0.402872 |

The baseline model predicts a larger total directional susceptibility for La in both fragments, but a **slightly larger increase for Dy when the second water is added**. These are different observables and are not contradictory. The same-origin zero-field dipole increments are approximately 0.617946 ea0 for La and 0.603011 ea0 for Dy; their ordering does not determine the susceptibility-increment ordering.

This conditional increment includes the added water's intrinsic and pair response and changes to the existing fragment. It is not alpha3, an isolated-water response, or a unique decomposition into induction/charge transfer. Only uniform lab-z response is measured. No agreement or disagreement with PBE0 is claimed until the matched quantum panel passes its own admission and field-step gates. The original full coordination-energy and zero-field trimer failures remain unchanged.

## Cost and preservation

Actual application: 9.662 s wall on one worker, 0.633 parent CPU-s and 7.297 aggregate native user+system CPU-s. Scheduler allocation accounting is root-owned. This saved review took 0.796 s wall/0.576 CPU-s on one local thread and made zero native calls.

A preliminary display attempted a nonexistent residual dictionary key after successful admission. Only that optional reporting extraction was corrected to the existing `maximum_printed_RMS_residual_Debye`; this is recorded in the review receipt. Native execution and all original values were untouched.
