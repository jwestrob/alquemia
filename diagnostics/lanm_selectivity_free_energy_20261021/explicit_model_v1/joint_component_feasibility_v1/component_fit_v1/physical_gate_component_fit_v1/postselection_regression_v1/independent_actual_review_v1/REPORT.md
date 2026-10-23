# New terminal consumed-trimer regression

**Saved arithmetic/native-readback review passed.** Six new trimer rows plus the admitted selected dependencies reproduce all42 raw terms, eight seven-term E3/mu3 cycles and four Dy-minus-La contrasts. Maximum arithmetic discrepancies are1.15e−14 kcal/mol and6.11e−16 Debye. No native, factor or solver call was repeated.

## Paired energy errors

Signed model minus retained electronic reference, kcal/mol. These are descriptive errors under the original PBE0/electronic-state limitations.

| Contrast | Baseline | Prior joint fit | Old component | New terminal |
|---|---:|---:|---:|---:|
| W8_near | -1.407051 | -5.311887 | -4.910487 | -3.558469 |
| W8_far | -1.039251 | -4.697730 | -4.312522 | -2.551999 |
| W9_near | -2.204283 | -5.771301 | -5.351390 | -4.115470 |
| W9_far | -0.878947 | -4.550055 | -4.164742 | -2.413251 |

All four paired energy errors improve from the old component terminal, but remain worse than baseline. Paired dipole-vector error norms (new; baseline in parentheses), Debye: W8_near 0.139094 (0.142977), W8_far 0.073223 (0.091906), W9_near 0.147676 (0.160159), W9_far 0.070345 (0.087935). All four improve versus every prior stage.

## Eight individual cycle tradeoffs

| Cycle | Old component E3 error | New E3 error | Old mu3 error norm | New mu3 error norm |
|---|---:|---:|---:|---:|
| W8 La near | -5.004212 | -5.083471 | 0.395804 | 0.400734 |
| W8 La far | -8.174940 | -8.396468 | 0.210559 | 0.216138 |
| W8 Dy near | -9.914699 | -8.641940 | 0.568231 | 0.539828 |
| W8 Dy far | -12.487462 | -10.948467 | 0.303943 | 0.289361 |
| W9 La near | -4.383229 | -4.508858 | 0.409095 | 0.416203 |
| W9 La far | -8.778273 | -8.995083 | 0.251671 | 0.257095 |
| W9 Dy near | -9.734619 | -8.624328 | 0.591752 | 0.563879 |
| W9 Dy far | -12.943015 | -11.408334 | 0.342712 | 0.327439 |

Energy columns are kcal/mol; dipole columns are Euclidean vector-error norms in Debye. All four Dy cycles improve from the old component terminal; all four La cycles worsen. Seven of eight individual absolute E3 errors remain worse than baseline (W9 La near is the exception). All eight individual dipole error norms remain better than baseline. The complete four-stage vectors, original print bounds and signed comparisons are retained in the machine-readable receipt and raw report.

## Interpretation and scope

The paired mechanical improvement survives in these consumed cycles, alongside improved paired dipoles. It does not remove the nonadditive energy discrepancy against the retained electronic reference. This is evidence for a remaining tradeoff; it is not proof of missing Hamiltonian terms, a global capacity limit, or a physical accuracy failure against a newly qualified reference. The eight cycles share objective dependencies and were consumed previously. Original reference and printed-moment uncertainty limits remain. No E3/mu3 acceptance threshold, physical admission or fresh-validation claim was introduced.

## Cost and evidence

Regression:15.8158 wall seconds;12.6223 parent CPU-s plus2.58151 native-child CPU-s, six new calls. Historical costs are not added. Independent review process:3.26954 wall /3.09948 child CPU-s, one thread, one invocation, zero native calls.

- [Review](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/postselection_regression_v1/independent_actual_review_v1/actual_v1/REVIEW.json) — `0327b6b83a44c2c63bfdc30c6dc79b75d55ccf5d343eb43edca07503d88a6021`.
- [Full four-stage report](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/postselection_regression_v1/independent_actual_review_v1/actual_v1/REPORT.txt) — `8712640f157a237496e6372b4a1f4baf4994c88936740acc0b1b97d1ef443cd7`.
- [Execution receipt](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/postselection_regression_v1/independent_actual_review_v1/EXECUTION_RESULT.json) — `5944e85e639f3996d424d2861daf7ec1e8add6cf33c350479ca04764688b7a01`.
- [Root source adoption](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/component_fit_v1/physical_gate_component_fit_v1/postselection_regression_v1/independent_actual_review_v1/ROOT_SOURCE_ADOPTION.json) — `7f88a74f9e203719f5836c64f844a6dce58951694c0c6979e7acd725705e02c3`.
