# Preliminary interpretation: more component flexibility, explicit response costs

**Snapshot only: five accepted steps, campaign still `started` when read.** The saved producer reports all 25 physical gates passing at these five incumbents. This is not terminal native admission. Final result/review and the new consumed-trimer result/review remain null in [VALUES.json](VALUES.json). No unfinished probe is interpreted, and no protocol, candidate or threshold is changed.

The present evidence supports a practical distinction: retaining the old fitted scores constrained the available component improvement, while relaxing those score ceilings exposes an electrical-response and energy-margin tradeoff. It does not establish global inadequacy of this Hamiltonian or physical acceptability of the new response errors.

## Absolute magnitudes to retain in the terminal report

P is the sum of literal metal–water interactions; B=C−P is the collective remainder of the metal–shell interaction. The RMS errors below concern the two **Dy-minus-La W8/W9 contrasts**, in kcal/mol, against the same consumed common-basis electronic-plus-original-D4 references. They are neither absolute parent energies nor exchange free energies.

| Quantity | Common origin24 | Previous constrained terminal | New accepted step5, preliminary |
|---|---:|---:|---:|
| P energy RMS error, kcal/mol | 26.6824 | 22.8146 | 15.0491 |
| B energy RMS error, kcal/mol | 22.3032 | 19.0105 | 12.3604 |
| Fixed-denominator minimax q | 1.000000 | 0.731099 | 0.318105 |
| Largest absolute original energy-control error, kcal/mol | 0.146086 | 0.167430 | 0.969274 |
| Largest rigid-body force-error RMS, kcal/mol/Å | 1.989753 | 1.816154 | 1.979287 |
| Largest projected body-component error, declared kcal/mol/Å | 4.910252 | 4.575280 | 4.345297 |
| Pair Cartesian dipole RMS error, D | 0.614949 | 0.614948 | 0.634584 |
| Existing C-dipole RMS error, D | 0.072309 | 0.072300 | 0.083650 |
| Susceptibility RMS error, a₀³ | 0.034161 | 0.031994 | 0.448326 |
| Original physical gates, producer-reported | 25/25 | 25/25 | 25/25 |

The force quantities are the existing weighted rigid-body projections, including their established torque convention; they are not atomwise Cartesian-force RMS. Their unchanged limits are RMS≤2 and maximum component≤5. Each of the ten energy controls retains |error|≤1 kcal/mol. At step5 the limiting energy row is balanced W9 exchange (−0.969274 kcal/mol); the largest force RMS is Dy/W9 (1.979287). Passing with small remaining margin is not zero error.

Electrical RMS values follow the existing context-balanced objective. No experimental acceptance threshold or uncertainty is supplied for these quantities. The susceptibility RMS increase is large in absolute terms as well as relative to the very small origin error; report 0.0342→0.4483 a₀³ rather than only a score ratio. The old E/F aggregate score likewise rises 0.153403→0.391129 despite all individual gates still passing.

Both component errors remain substantial. At step5 the signed errors are P8=+14.1241, B8=−10.7183, P9=+15.9204, B9=−13.8086 kcal/mol. Opposite signs still cancel in C=P+B, giving derived common-basis C errors +3.4058/+2.1118. These derived sums are not extra observations or substitutes for the original physical-gate reference cycles. Relative to the original pre-water-fit component baseline (P/B RMS 13.1700/13.0947), B is now smaller but P remains larger. The declared combined recovery milestone and ≤1 kcal/mol component development target are both still unmet at this snapshot.

## What separates preservation policy from physical tradeoff?

The clean policy comparison is the already admitted **matched local conic comparison**: same origin, Jacobian, trust region, parameter box and physical inequalities, removing exactly four old-score cones. Its full strict witness improved predicted q from the full-constraint witness’s 0.842764 to 0.438713. Solver status stayed `AlmostSolved`; feasibility did not relabel convergence. That establishes available local directions outside the preservation ceilings, not an optimal global solution.

The single full-fraction native pilot realized q=0.442386 and P/B RMS 17.3193/14.8343, but failed the water-minus energy gate at −1.002200211 kcal/mol. That failed point remains excluded. The current finite campaign’s accepted smaller steps demonstrate that useful component improvement can coexist with the unchanged physical gates, while the response scores worsen. The complete terminal audit must confirm this before admission.

Comparing the previous eight-step SLSQP campaign with this Clarabel campaign alone does **not** isolate a pure ceiling effect: solver formulation and subsequent trajectories also differ. Use the matched local comparison to identify the policy effect and the nonlinear trajectories to show its practical consequences. Do not attribute the entire difference between terminals to one removed block, or claim that a finite stop proves no better point exists.

## Consumed-trimer comparison awaiting the new terminal

The planned regression retains all eight seven-term E3/μ3 cycles (W8/W9 × La/Dy × near/far) and all four paired contrasts. It reuses immutable reference/baseline data and exact model states; at most six missing analyze programs are allowed. These previously consumed cycles share dependencies with development data. They test behavior beyond the selected P/B objective, **not independent holdout transfer**.

The previous constrained terminal improved all four paired E3 discrepancies over the preceding joint-water fit, but all four remained worse than the original pre-water-fit baseline:

| Dy−La E3 signed error, kcal/mol | Original baseline | Previous constrained terminal | New terminal |
|---|---:|---:|---|
| W8 near | −1.407051 | −4.910487 | pending |
| W8 far | −1.039251 | −4.312522 | pending |
| W9 near | −2.204283 | −5.351390 | pending |
| W9 far | −0.878947 | −4.164742 | pending |

Across its eight individual cycles, previous-terminal E3 errors ranged from −12.9430 to −4.3832 kcal/mol; seven were worse in absolute magnitude than the original baseline. Individual μ3 vector-error norms were 0.21056–0.59175 D, all smaller than baseline, while all four paired μ3 errors remained worse than baseline (0.09104–0.18266 D). Energy and dipole behavior therefore need separate reporting. [VALUES.json](VALUES.json) retains all eight comparator rows and four paired rows; no favorable subset is chosen.

## Terminal decision structure — descriptive, no new gate

1. **Establish the actual selected state.** Bind the frozen terminal result/configuration/selected artifact to the passed terminal-review receipt. Preserve stopping reason, literal solver statuses, all failed or untouched trials, selected46/native provenance, actual counters and measured costs. A failed or unavailable audit supplies no candidate for the regression.
2. **Report component repair and its price.** Replace the preliminary column with terminal P8/B8/P9/B9, P/B RMS, derived C sums, q, all 25 physical gates and their margins, four old scores, and absolute dipole/susceptibility RMS. Show origin24, previous constrained terminal and original pre-water-fit component baseline distinctly. Lower q alone is not small absolute error or electrical validation.
3. **Report all eight regression cycles and four contrasts.** Bind its actual review; retain model/reference values, signed and absolute E3 errors, μ3 error vectors/norms, print bounds, reference flags and reused/new-call counts. Compare with both previous terminal and original baseline without introducing E3/μ3 acceptance thresholds.
4. **Interpret the observed combination.** P/B improvement with surviving physical gates supports flexibility beyond the old ceilings. If trimer discrepancies also improve, report the extent and any remaining baseline deficit; if they worsen, report a finite aggregate-versus-three-body tradeoff. If electrical errors grow, retain that cost even when trimer energies improve. None of these cases identifies a unique missing interaction term or establishes a free-energy preference.

The reference-method, frozen electronic-state, finite geometry and consumed-data limitations remain. B includes all nonpair orders and is not a unique polarization/charge-transfer observable. No fresh radial, donor/protein, periodic dynamics or held-out outcome is opened or released by this synthesis.

## Evidence and cost

[VALUES.json](VALUES.json), SHA256 `ca9a8b476fa608fe9c114fe85f0e6e7ac65907210075400861fa2fe640542bfd`, records exact accepted-artifact pins, the live-result hash at one read, prior-terminal/result/review authorities, source declarations and the deliberately null final bindings. Extraction cost: 0.205 s wall / 0.167 CPU-s, one thread, zero model/reference/optimizer calls. This interpretation does not repeat native admission or amend shared checkpoints.
