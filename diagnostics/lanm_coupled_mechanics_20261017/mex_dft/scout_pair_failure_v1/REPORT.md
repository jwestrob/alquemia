# Mex embedded-DFT scout fails directional force consistency

17 October 2026 UTC. The original La external-charge displacement pair is complete. Saved collector 13271 replays both native endpoints and returns `failed_derivative_gate`. This is a completed scientific comparison with a failed acceptance gate; the scheduler exit 1 is not an infrastructure crash.

## Result

Derivatives along the frozen A/97//O displacement, in kcal/mol/Å:

| Component | Analytic dE/ds | Centered energy difference | Difference | Allowed difference | Pass |
|---|---:|---:|---:|---:|---|
| Embedded PBE0 | -3.673754 | -4.087955 | -0.414201 | 0.005409 | No |
| Retained MM | 17.752465 | 17.752485 | 0.000020 | 0.006775 | Yes |
| Fixed-charge ALPB | -5.530643 | -5.530642 | 0.000001 | 0.005553 | Yes |
| Total | 8.548069 | 8.133889 | -0.414180 | 0.005855 | No |

The discrepancy is localized to the embedded PBE0 contribution in this directional test. It prevents qualification of the assembled force under the original criterion. It does not establish a sign/mapping defect, an inaccurate physical DFT model, or which alternative model is accurate.

Both singlepoints terminate normally and are available under the unchanged native admission policy. The origin and displaced core coordinates are byte-identical; the declared external-charge position changes. The original template and endpoint template differ only by the declared removal of EnGrad. Source/state identities and the finite-difference denominator pass the original collector. All eight declared endpoint slots remain present; the other six were not requested.

## Numerical and electronic limitations

Origin, plus and minus terminate via energy, gradient and energy checks, respectively. All retain unmet printed residual criteria and unqualified physical states. Native convergence is not a proof of electronic stability or derivative accuracy. The original origin used 112 MPI ranks/MaxCore10800; current endpoints use 12 MPI ranks/allocation-derived MaxCore3100. This runtime difference is preserved, not established as a cause.

The collector subtracts exact printed DFT energy strings in native Hartree. Its printed-energy and binary64 sensitivity scales are below 0.000001 kcal/mol/Å. The last-SCF-energy-change proxy is 0.001920 kcal/mol/Å, but it is not an energy-error bound. The fixed 0.001 Å test has no step-halving or truncation qualification. These observations do not identify the cause of the 0.414201 discrepancy.

## Decision and next work

The original failure rule stops the six remaining qualification singlepoints and the dependent eight-cell Mex response panel. No changed solver, repeat, alternate step size or favorable subset is launched. Saved source, field-gradient mapping and native-state evidence are being reviewed to identify an actionable mechanics issue. The independent whole-peptide before-La force reference 13051 continues; it has no external point-charge field and is a separate test.

## Actual costs and preserved evidence

| Job | Work | Elapsed seconds | CPUs | Allocated CPU-seconds |
|---|---|---:|---:|---:|
| 12977 | Plus endpoint | 26,666 | 12 | 319,992 |
| 12978 | Minus endpoint | 22,132 | 12 | 265,584 |
| 13271 | Saved pair collection | 6 | 1 | 6 |
| Total | This scout and collection | | | 585,582 |

The two endpoints make two DFT singlepoints, two MM configurations/eight group queries and two ALPB calls. The collector makes zero molecular calls. These costs exclude previously completed preparation and reused origin chemistry. Plus final gridx3 took 1,875.779 s; minus took 1,920.673 s.

The exact collector config, declaration, submission, result, stdout/stderr and scheduler accounting are archived here. ROOT_REVIEW.json pins the authoritative workspace result and native endpoints; NATIVE_REVIEW.json records independent review. The plus raw archive is ../scout_plus_result_v1/; the earlier minus archive remains ../scout_minus_result_v1/. Original failures, outputs and gates are unchanged. No preference classifier, production change or completed-goal claim follows.
