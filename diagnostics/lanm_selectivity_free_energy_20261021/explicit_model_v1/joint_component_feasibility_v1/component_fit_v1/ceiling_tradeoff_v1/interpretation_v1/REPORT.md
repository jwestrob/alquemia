# Preservation ceilings materially restrict the local response

Removing only the four development score ceilings unlocks a much better **linearized** P/B compromise at the unchanged origin, trust region and physical thresholds. The full returned physical direction has actual minimax q **0.4387129496723387**, versus the certified full-constraint baseline **0.8427637923620422**: an absolute gain of **0.4040508426897035**. This is strong evidence of a local preservation-policy tradeoff. The restricted campaign alone cannot support a conclusion that the Hamiltonian lacks response flexibility.

The result remains **AlmostSolved**, after 17 iterations. Its gap in t is 1.22962e−10, just above the frozen 1e−10 requirement; scaled primal/dual residuals are 7.03167e−13 / 6.98631e−15. No status, tolerance or gate is changed. Witness feasibility is separately supported by the original scalar physical checks, with only auxiliary epigraph bookkeeping. No exact nonlinear candidate has been evaluated.

| Existing fraction | Actual q | Gain against baseline | Strict retained feasibility | Frozen ≥0.05 trigger |
|---|---:|---:|---|---|
| 1 | 0.438713 | +0.404051 | Pass | Pass, first |
| 1/2 | 0.690633 | +0.152131 | Pass | Pass |
| 1/4 | 0.838136 | +0.004628 | Pass | Fail |
| 1/8 | 0.917273 | −0.074509 | Pass | Fail |

All25 physical gates pass in these predictions. The full fraction sits very close to boundaries: aquo9 balanced-exchange error is 0.9999999999985 kcal/mol against 1, paired water-plus force RMS is 1.9999999999974 against 2 kcal/mol/Å, and minimum box margin is 9.86e−14. The minimum retained scalar margin is 1.48526e−12. These are strict ordinary-double certificates with little physical headroom; they do not establish nonlinear feasibility. The raw full-point epigraph already passes by 5.02e−14; the declared auxiliary cushion is retained uniformly and excluded from every gain calculation.

## The price in previously fitted scores

All four old score ceilings are exceeded by the full fraction. These are objective scores, not four newly failed experimental tolerances.

| Old block | Origin ceiling | Predicted full-fraction score | Ratio to ceiling |
|---|---:|---:|---:|
| E/F | 0.1534031 | 0.4690138 | 3.05739 |
| Pair dipole | 0.3781620 | 0.3838325 | 1.01499 |
| Collective dipole | 0.00522866 | 0.00539378 | 1.03158 |
| Susceptibility | 0.00116701 | 0.0440068 | 37.7091 |

The P/B RMS predictions are approximately 17.25 / 14.77 kcal/mol; substantial component errors remain. The large relative susceptibility loss and aggregate E/F loss must be reported even if all original physical gates survive. Small dipole-score increases are also retained. No objective improvement makes those tradeoffs physically acceptable by itself; no independent μ/χ accuracy tolerance is invented here.

## Frozen policy evaluation — report only

Policy `NONLINEAR_FOLLOWUP_POLICY.json` SHA256 `9e9949519d32e935511016a1c2d4a484e673ba0c90081f700e607c3c6a65a6fa` selects the first strictly feasible existing fraction with gain ≥0.05. It therefore identifies **fraction 1**, key **d1fc95754c0cd8dceb9a7616caf2e8409974f6952d0a904e0f2ead9645acbec3**. This is a prospective development trigger, not a physical qualification or execution release. The half fraction is not substituted after observing the full fraction's boundary proximity.

The useful next check is the separately declared, single-vector nonlinear evaluation specified by that policy, contingent on independent actual arithmetic review, source preparation and root release. It should test whether the large predicted gain survives while reporting all25 physical gates and all four old-score losses. A native failure or gate failure must remain a failure; no second fraction, optimizer or automatic fit follows. No native action is taken by this report.

## Evidence and cost

Saved producer RESULT SHA256 `1d4c5db9a34a5baf2ad93dd9168b30432c151091a8b46f0cc1b8798247fb38d6`; [VALUES.json](VALUES.json) pins result, configuration, source review, release, baseline and policy and preserves all four exact parameter vectors/predicted scores. Independent actual arithmetic review was pending when this interpretation was written. Producer cost: **1.56752 s wall / 1.49989 CPU-s**, one solver call, one numerical thread, zero native/reference/candidate evaluations. This report invokes no solver or native program.
