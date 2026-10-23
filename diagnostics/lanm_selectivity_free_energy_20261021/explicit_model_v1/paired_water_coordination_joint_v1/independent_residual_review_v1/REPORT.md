# Independent interpretation of terminal joint-water fit

**The evidence supports local least-squares stationarity with substantial residual error. A simple iteration cap, parameter bound or trust-step restriction does not explain this stop.** The original scientific failure remains: 9/10 energy and 13/15 force controls pass. Independent raw review17683 replays all ten selected native/reference cells and the exact objective, while retaining scientific acceptance as false.

## Quantitative evidence

- Job17676 stops after a complete non-improving batch at step22 of24, using338 of384 new vectors and4120 of5120 native programs.
- All six parameters are interior. The nearest normalized bound is0.04237 away; the proposed unconstrained step is at most2.57×10⁻⁶, versus the0.025 trust limit. No final parameter probe is clipped.
- The last Jacobian origin equals the selected point. An independent reconstruction from all twelve saved derivative vectors reproduces the413×6 matrix exactly. SVD at the frozen cutoff retains all six columns; condition number149.169.
- The unrestricted local linear projection reduces weighted loss by only1.55×10⁻¹⁰, from0.895803723689 to0.895803723534. All four evaluated final trials instead increase it by1.4–2.8×10⁻⁹; none passes the development gates.
- The W9 exchange error remains−4.338284kcal/mol; its projected value is−4.338308. This one residual contributes65.66% of the aggregate weighted loss. The other failures are La/W9 rigid RMS2.042854 and paired/dimer-plus RMS2.203824, each exceeding the unchanged2 threshold.

## Interpretation and limits

The available six-parameter linear response cannot appreciably reduce the current weighted residual. Simply continuing the same local procedure, increasing its iteration cap, or loosening bounds has no support from these saved results.

This is **not a feasibility certificate for all individual thresholds**. A weighted least-squares stationary point can fail individual gates even if another, potentially higher-loss point satisfies them. The analysis also does not exclude distant nonlinear solutions or establish a global family limitation. The Jacobian uses finite normalized halfsteps0.005; tiny final loss differences may reflect finite-difference or printed-output precision. The multi-kcal/mol W9 discrepancy remains much larger than those final changes.

The residuals are correlated engineering controls, not independent observations defining statistical confidence. The projection identifies a limitation of the local parameter response; it does not identify a unique missing physical interaction. W8 and W9 also differ in shell geometry, so the failure cannot be attributed solely to coordination number.

The source-only Cartesian and fresh radial adapters remain unbound and held by the failed scientific gate. No parameter update, extension, new molecular call or fresh validation read followed this review.

## Evidence

The matching workspace `independent_residual_review_v1/RESULT.json` (SHA256 `c8c0535469e8e49c013307b3239762d7a4d3f300d84dcf370332dc553864acdc`) records23 passing independent checks, exact parent pins and the reconstructed projection. Saved arithmetic used one thread,2.935s wall/2.931CPU s, zero native calls, and did not repeat objective or raw-native replay.
