# Fixed direct/mutual response sensitivity: physical interpretation

All four profiles executed and were independently admitted. **None resolves the consumed panel.** Separate direct and mutual damping produce distinguishable responses, but the observed improvements are small relative to the collective discrepancies and have clear tradeoffs. This is a finite sensitivity result, without a fitted or selected candidate.

## Evidence and scope

Job 18413 completed 120 cells with 236 new native programs, four exact saved gradient reuses, zero new Cartesian finite differences, and no application error. Application wall time was 44.339 s. Root's terminal scheduler record reports 48 s on four CPUs, 192 allocated CPU-s, TotalCPU 145.272 s, and MaxRSS 1,493,344 K.

The frozen profiles apply the same increment to both metal–water-O pairs: D± changes only gamma_D by ±5%; M± changes only gamma_M by ±5%. Geometry, radii, epsilon, polarizabilities, other pair interactions and the original thresholds remain fixed. Increasing gamma makes the cubic damping factors less screened at a fixed separation; it does not simply rescale every molecular response uniformly.

Independent actual review `d3c624ab02e1931905c6cb53ce7ccc582df0bf13187c0b20726a0d93bd81c379` admitted all raw rows, receipts and summary arithmetic. The saved reviewer normalized only four proven tuple/list representations of profile increments; the executed source and result remain unchanged. Maximum admitted induction residual was 9.1e-9 and analyze/gradient energy closure was exact at printed precision.

The accompanying workspace `physical_review_v1/SUMMARY.json` retains all 25 numeric controls for the baseline and every profile, all 32 trimer energy/dipole outcomes, all 16 center outcomes, and exact input/result/review pins. Its SHA256 is `450311d99f0d1a7f093b45926ccbd6383b88a3814e90a36c36776546731381b9`. The application result SHA256 is `9179aa66786acc311de01bc9f98e3346a1b29f9e6eb6011d78828509216982e8`. This interpretation used saved data only and made zero native calls.

## Original 25 checks

Energy gates retain |error| ≤1 kcal/mol. Each force fixture must satisfy rigid RMS ≤2 and maximum body component ≤5 kcal/mol/Å. Pass counts include every original slot.

| Profile | Objective | Energy /10 | Force /15 | Total /25 | Original W9−W8 paired error, kcal/mol |
|---|---:|---:|---:|---:|---:|
| Baseline | 0.895804 | 9 | 13 | 22 | −4.338284 |
| D+ | 1.122275 | 7 | 13 | 20 | −4.349483 |
| D− | 1.116930 | 9 | 9 | 18 | −4.321582 |
| M+ | 0.888390 | 9 | 13 | 22 | −4.328993 |
| M− | 0.903286 | 9 | 13 | 22 | −4.347440 |

- **D+** gains the La/W9 force pass: RMS 2.042854→1.870026. It loses the compressed La–water force pass (RMS 1.014845→2.418871) and the compressed/central single-water exchange-energy passes (errors −1.180387/−1.170146). The paired expanded-water force failure remains. Objective rises 25.28%.
- **D−** gains no pass and loses four force passes: La/expanded water, Dy/compressed water, Dy/W9, and paired/compressed water. Their RMS errors become 2.641992, 2.375883, 2.175117 and 2.019761; La/expanded water also exceeds the maximum-component gate (5.236035). Objective rises 24.68%.
- **M+** changes no pass/fail outcome. Objective decreases 0.83%; the W9 energy error improves by only 0.009290 kcal/mol. The original failures—W9 paired energy, La/W9 force, paired expanded-water force—remain.
- **M−** also changes no gate and raises the objective 0.84%. Its small local force improvements do not compensate for the worsening W9 contrast and collective moment response.

All six individual single-water force fixtures still pass under M±, although their numerical errors shift in both directions. This is useful evidence that the mutual response can be changed without a new isolated-water force failure at these amplitudes. It does not establish a successful repair: even the largest favorable change in the original coordination error, D−'s +0.016702 kcal/mol, removes only 0.39% of its magnitude.

## All eight trimer increments

E3 is the unchanged seven-term M+Wi+Wj nonadditive energy. Reference values here are the existing **common-parent-basis electronic** quantities; they are distinct from the original total-energy development targets above. No target was replaced.

| Profile | E3 errors improved /8 | Largest absolute E3 change, kcal/mol | Remaining E3 error range, kcal/mol | Dipole vector errors improved /8 | Largest magnitude change, Debye |
|---|---:|---:|---:|---:|---:|
| D+ | 8 | 0.122282 | −8.129301 to −4.504697 | 8 | 0.002748 |
| D− | 0 | 0.134164 | −8.246294 to −4.743521 | 0 | 0.003091 |
| M+ | 4 | 0.084913 | −8.144340 to −4.638206 | 8 | 0.002809 |
| M− | 4 | 0.087828 | −8.226380 to −4.594639 | 0 | 0.002915 |

D+ raises E3 in every context by 0.054–0.122 kcal/mol and strengthens every collective dipole magnitude slightly. That helpful absolute trend coexists with worse paired metal discrimination: La changes more than Dy. D− reverses those absolute trends, despite its small improvement in the full-cluster coordination contrast.

The mutual response has an angle-dependent tradeoff: M+ improves all four **far** E3 errors but worsens all four **near** errors; M− does the reverse. This difference from D's same-sign E3 response establishes distinct local response directions in the tested panel. It does not justify adding the two finite perturbations or extrapolating a coefficient beyond the tested amplitudes.

Every profile still underestimates every collective dipole magnitude. Across this entire fixed matrix the shortfalls remain 58.47–85.92%. The paired magnitude ordering remains wrong in every one of the four contexts:

| Context | Reference | Baseline | D+ | D− | M+ | M− |
|---|---:|---:|---:|---:|---:|---:|
| W8 near | +0.120781 | −0.022197 | −0.023085 | −0.021235 | −0.023324 | −0.020992 |
| W8 far | +0.064076 | −0.027830 | −0.028822 | −0.026747 | −0.029482 | −0.026173 |
| W9 near | +0.139398 | −0.020761 | −0.021853 | −0.019607 | −0.022054 | −0.019380 |
| W9 far | +0.060279 | −0.027656 | −0.028617 | −0.026598 | −0.029306 | −0.025998 |

Entries are |mu3(Dy)|−|mu3(La)| in Debye. Positive means stronger Dy suppression. The reference has positive values throughout, while every model value remains negative. Positive D/M increments modestly improve individual magnitudes but make this ordering discrepancy larger; negative increments slightly reduce the ordering discrepancy while weakening already insufficient individual magnitudes. Both far E3 metal orderings also remain reversed under all four profiles. Near E3 ordering retains the reference sign but remains too weak.

## Four cluster centers and interpretation

W9 interaction C-dipoles (Debye, z component; transverse components are zero for the model):

| Metal | Reference | Baseline | D+ | D− | M+ | M− |
|---|---:|---:|---:|---:|---:|---:|
| La | 2.244444 | 2.107200 | 2.126975 | 2.085079 | 2.113334 | 2.100869 |
| Dy | 2.139823 | 2.017552 | 2.026776 | 2.006823 | 2.020143 | 2.014752 |

Positive increments help the already modest full-cluster dipole underresponse, again more for La. W8 total model C-dipoles remain zero by symmetry under every profile; that cancellation is not evidence of physical accuracy. M± leave the W8 center forces unchanged at printed precision, consistent with the baseline's symmetry-suppressed metal induced dipole, while its trimer and W9 responses remain active. D± change the W8 forces as well.

**Mechanical conclusion:** the separate direct/mutual controls are numerically active and distinguishable. D+ gives a small, consistently favorable absolute trimer response at the cost of isolated-force/energy agreement and worse metal ordering. M+ gives a gentler change with no new original gate failures, but its near/far energy tradeoff and unchanged wrong metal ordering leave the central collective-response problem unresolved. The present evidence supports neither a repaired model nor another optimizer extension. It also does not show that every possible parameter value or functional form must fail.

Only the declared 12 extracted dimers were recomputed per profile, so this panel does **not** evaluate full-shell P or B under the modified kernels. The trimer energies/dipoles are specific nonadditive observables, not a unique decomposition into physical induction, charge transfer or penetration. No fitting, new quantum calculation, held-out evaluation, model promotion or proposed submission is part of this review.
