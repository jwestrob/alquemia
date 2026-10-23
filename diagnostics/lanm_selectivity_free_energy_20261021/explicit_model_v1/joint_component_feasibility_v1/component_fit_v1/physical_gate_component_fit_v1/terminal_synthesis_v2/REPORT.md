# Terminal synthesis: component recovery with unresolved electrical and three-body errors

The unchanged twelve-parameter Hamiltonian reaches **P/B RMS errors of 12.7203/10.4418 kcal/mol and passes all 25 original energy/force gates**. Both component errors improve over the original pre-water-fit model. This establishes useful flexibility beyond the former four score-preservation ceilings. It does **not** establish adequate absolute accuracy, electrical transfer, or readiness for free-energy sampling.

The selected point is `feca1810…`, admitted by fit review `13f23ffd` against terminal result `4a374ff3`. The independent consumed-trimer review `0327b6b8` admits arithmetic/native provenance for result `686d5716`; it introduces no physical-accuracy threshold. [INPUTS.json](INPUTS.json) pins these authorities, the selected artifact, frozen objective/domain and preserved preliminary reports. [VALUES.json](VALUES.json) retains every quantity summarized below.

## Component recovery and remaining cancellation

P is the sum of literal metal–water interaction cycles; B=C−P is the collective remainder of the metal–shell cycle. These errors concern the **two Dy-minus-La W8/W9 contrasts**, against the consumed common-basis electronic-plus-original-D4 references. B includes all nonpair orders; it is not uniquely a polarization or charge-transfer contribution.

| Component RMS error, kcal/mol | Original pre-water-fit model | Origin24 | Previous ceiling-constrained terminal | New terminal |
|---|---:|---:|---:|---:|
| P | 13.16996 | 26.68241 | 22.8146 | **12.72032** |
| B | 13.09468 | 22.30315 | 19.0105 | **10.44183** |

Relative to origin24, P/B RMS decreases 52.33%/53.18%; relative to the original model, the gains are **3.41%/20.26%**. Fixed-denominator minimax q reaches 0.22727182. The prospectively declared recovery milestone passes, while the ≤1 kcal/mol component development target remains unmet. Neither the normalization nor that development target is a reference uncertainty estimate.

| Signed paired error, kcal/mol | W8 | W9 |
|---|---:|---:|
| P | +11.93611 | +13.45891 |
| B | −8.92405 | −11.76541 |
| Derived C=P+B | +3.01206 | +1.69350 |

Opposite component errors still cancel. The smaller total C error therefore cannot stand in for accurate pair and collective interactions. These common-basis sums are distinct from the original balanced-exchange physical-gate cycles. Missing per-parent D4-complete targets are not invented.

## Electrical errors: retain signs and conditional differences

| Existing objective quantity | Origin24 RMS error | Terminal RMS error | Consumed reference RMS scale |
|---|---:|---:|---:|
| Pair Cartesian dipoles, D | 0.614949 | **0.639432** | 2.271140 |
| Full-shell C Cartesian dipoles, D | 0.072309 | **0.085846** | 0.895190 |
| Four susceptibility states, bohr³ | 0.034161 | **0.680444** | 17.975770 |

These are the original context-balanced weights and absolute units: four equally weighted parent contexts for pair/C, and four equally weighted susceptibility states. The block scores precede the outer factor 1/4. Scales are engineering choices, not physical tolerances. The old E/F score rises 0.153403→0.386590; all four old score ceilings are exceeded descriptively, as permitted in this campaign.

Susceptibility uses the frozen reference energy-curvature estimator at h=0.001 atomic units:

| State | Reference, bohr³ | Model, bohr³ | Model−reference, bohr³ |
|---|---:|---:|---:|
| La dimer | 15.47282170 | 14.59369477 | **−0.87912693** |
| La trimer | 22.01652416 | 21.07746852 | **−0.93905564** |
| Dy dimer | 13.15810670 | 13.49032327 | **+0.33221658** |
| Dy trimer | 19.88069016 | 20.17557433 | **+0.29488416** |

The same-metal errors largely cancel in the **conditional trimer-minus-dimer increment**. La gives model/reference 6.48377375/6.54370245 bohr³, error −0.05992871; Dy gives 6.68525105/6.72258346, error −0.03733241. Their Dy-minus-La difference is **0.20147730 versus 0.17888101 bohr³**, error +0.02259629. Thus the aggregate state RMS does not imply equally poor incremental response. These derived differences are not new independent targets or accuracy admissions.

The full-shell C dipole has a distinct sign problem. W9 model z components are La 1.95297050 D and Dy 2.08087005 D, with signed errors −0.29147358 and −0.05895607 D. Consequently **model Dy−La is +0.12789955 D, while reference Dy−La is −0.10461795 D**: contrast error +0.23251751 D. W8 model z is effectively zero, versus reference 0.00096555/0.00083266 D; all four x/y error components are below 8.55×10⁻⁷ D in magnitude. The complete signed vectors are saved. Percentages of nearly cancelling W8 or transverse components would be misleading.

## Consumed trimer regression: energy and dipole outcomes differ

All eight seven-term E3/μ3 cycles and four paired contrasts are included. Six new analyze programs were needed; unchanged dependencies were reused. These are consumed regressions sharing development dependencies, not holdouts. Their retained **electronic E3** references must not be relabeled D4-complete total-energy targets.

| Dy−La E3 signed error, kcal/mol | Original baseline | Prior joint fit | Previous component terminal | New terminal |
|---|---:|---:|---:|---:|
| W8 near | −1.40705 | −5.31189 | −4.91049 | **−3.55847** |
| W8 far | −1.03925 | −4.69773 | −4.31252 | **−2.55200** |
| W9 near | −2.20428 | −5.77130 | −5.35139 | **−4.11547** |
| W9 far | −0.87895 | −4.55006 | −4.16474 | **−2.41325** |

All four paired energy errors improve over both fitted comparators but remain worse than the original baseline. Among eight individual cycles, all four Dy energy/dipole errors improve over the previous component terminal and all four La errors worsen. Seven of eight absolute E3 errors remain worse than baseline; terminal signed errors range from −11.40833 to −4.50886 kcal/mol. Conversely, all eight individual μ3 error norms improve over baseline, and all four paired μ3 error norms improve over every comparator (new 0.139094/0.073223/0.147676/0.070345 D in table order). [The complete eight-cycle interpretation](../postselection_regression_v1/independent_actual_review_v1/REPORT.md) preserves raw-token/printing limitations and every comparison.

## What the finite fit does and does not establish

The tightest actual margins remain small: balanced W9 exchange error −0.99384220 kcal/mol leaves **0.00615780** to the original ±1 gate; Dy/W9 rigid-body force RMS 1.99487845 leaves **0.00512155 kcal/mol/Å** to its limit 2; paired water-plus RMS leaves 0.00862151. The largest projected body-component error is La/W8 4.89593512, leaving 0.10406488 to limit 5. These are the established body-force/torque projections, not atomwise force RMS.

The parameter box is also relevant. **La_u=0.93498087** is 0.06501913 below its upper bound; **Dy_v=−0.99965881** is only 0.00034119 above its lower bound. These are dimensionless box coordinates, not native c/d. With c=S(u+v)/4 and d=S(u−v)/4, native La c/d=0.51370035/0.13437901 and Dy c/d=−0.11556035/0.57735034. Box limits were declared development bounds, not evidence of physical saturation.

All eight outer steps completed; the campaign stopped at its **step limit**, with q still decreasing 0.25507470→0.22727182 in the last step. Seven proposal records say `AlmostSolved`, one says `Solved`; independent witness checks and actual native gates governed acceptance without relabeling solver convergence. There were 226 complete vectors (199 new, 27 imported), 6,440 new native programs and seven new solver calls plus one imported proposal. The last Jacobian is at `a879be2f…`, the selected step7 point, **not final `feca1810…`**. This is a finite developmental result, with no global optimality, local final-point infeasibility, or Hamiltonian-capacity proof.

## Recommended next decision: one final-point local attribution panel

Before adding model terms or continuing calibration, distinguish domain/gate restrictions from electrical tradeoffs at the **exact final point**. First inventory the existing caches and use the saved step7-origin Jacobian/trials only as off-center evidence. Then, under a separate finite declaration, obtain only missing final-point ±0.005 normalized probes (clamped to the unchanged box; record their actual, potentially asymmetric spans). The ceiling is 24 vectors/1,344 programs before exact per-metal reuse, with no new references or trimer fitting.

Use one fixed linearization for a small, preregistered comparison of paired P/B improvement under the original 25 gates and parameter domain, then with explicit signed electrical-error nonincrease rows. Report which domain faces, individual physical gates and electrical rows restrict feasible local improvement; diagnostic relaxation of a family may identify sensitivity but cannot admit a point violating that family. Electrical nonincrease would be a clearly labeled development comparison, not an invented experimental tolerance. Independent primal margins and actual saved trial prediction errors must accompany solver status. No solver/native work is released by this note.

This panel can identify a specific constraint or residual direction worth addressing. A failed local affine comparison remains local and may reflect derivative/trust-region or numerical limitations. Successful local directions would motivate one separately declared nonlinear check. Neither outcome alone warrants new terms, a wider domain, more fit steps, or production sampling. Reference-method/frozen-state limitations remain, and periodic engineering qualification is a separate branch.

## Cost and evidence boundary

The literal terminal **`fit.cost`** records 2,012.397437 s wall and **5,086.966540 measured CPU-s** = 2,112.860558 controller + 1,267.453806 child user + 1,706.652176 child system. Execution was truthful login mode with no Slurm allocation; inherited costs are excluded. The root receipt's historical `fit_actual_cost:null` is not used. Trimer execution took 15.815828 s wall, 12.622259 controller CPU-s and 2.581510 child CPU-s. This synthesis extraction took 0.259129 s wall/0.246957 CPU-s on one thread, with zero native, solver, reference or heldout calls. It relies on existing admission reviews and changes no model, threshold, source, or checkpoint.
