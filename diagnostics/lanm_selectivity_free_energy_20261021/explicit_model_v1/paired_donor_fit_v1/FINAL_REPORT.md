# Shared donor fit: improved, but fails development requirements

21 October 2026. Final result of job 17621; independently replayed against saved native outputs.

**The fit reduces its weighted objective by 68%, but still fails 21 of 33 energy
checks and 27 of 36 force fixtures. This candidate is not qualified for protein
or aqueous exchange production.**

The question was whether eight shared metal–donor parameters could reproduce
the consumed La/Dy carboxylate and amide data while retaining the previously
passing water model. Mono- and bidentate acetate share a carboxylate parameter
pair for each metal. The water model and all reference targets stayed fixed.

## Actual outcome

| Check | Initial failures | Final failures | Largest final error |
|---|---:|---:|---:|
| Aquo-referenced exchange energy | 11/12 | 10/12 | 3.949 kcal/mol |
| Wet-minus-dry exchange energy | 3/3 | 3/3 | 5.805 kcal/mol |
| Individual radial work | 9/12 | 7/12 | 2.524 kcal/mol |
| Paired radial work, control only | 3/6 | 1/6 | 1.995 kcal/mol |
| La force fixtures | 12/12 | 11/12 | 15.585 kcal/mol/Å RMS |
| Dy force fixtures | 12/12 | 8/12 | 13.459 kcal/mol/Å RMS |
| Dy−La force fixtures | 12/12 | 8/12 | 9.402 kcal/mol/Å RMS |

The unchanged limits are 1 kcal/mol per energy check and, for each force fixture,
RMS ≤2 plus maximum component ≤5 kcal/mol/Å. A force fixture can fail its maximum
component even when its RMS passes. The objective changed from 28.6273422 to
9.14353359 over 20 steps, 376 new parameter vectors and 6,528 native program
invocations. The final complete trial batch produced no actual improvement.

## What this resolves

The numerical improvement is real: independent raw-input, native-output and
objective replay passes 247 checks. All eight parameters are inside their bounds.
The final measured response matrix has full rank 8/8 and condition 72.87.
At this same selected point, its unrestricted local linear estimate lowers the
objective by only about 0.000000018; actual half and full trial steps do not
improve it. More iterations of this same local fit are not supported.

The remaining errors include the change caused by adding surrounding waters and
large donor-force errors. They motivate testing the interaction response itself.
They do not identify a unique missing term or prove that every distant parameter
choice in this model family must fail.

These are consumed development configurations and provisional frozen-core
ECP-PBE0-D4 targets. Reference uncertainty remains. No reserved donor validation
was opened, no new fit was started, and no exchange free energy is established.

## Evidence and preserved reporting failure

- [Independent actual review](../../../../workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/actual_review_v1/result_v1/REPORT.md) lists every energy and force check.
- [Local response analysis](../../../../workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/span_diagnostic_v1/result_v1/REPORT.md) records the predicted and actual final trials.
- [Figure, PNG](../../../../workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/report_recovery_v2/report/donor_response.png) and [PDF](../../../../workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/report_recovery_v2/report/donor_response.pdf) show the errors, forces and optimization history. The PNG was visually inspected; labels and panels are legible.

Report job 17625 wrote its text, then failed because its Python interpreter lacked
Matplotlib. The failed attempt is retained. The unchanged report source and result
were rendered into a new directory with the existing Python 3.12 environment;
the recovery used one local thread, took 22.027 seconds and made zero molecular
calls. This was an export repair, not a repeated fit.

Scheduler cost: fit 143,520 allocated CPU-seconds; failed export 11; actual review
17; local-response diagnostic 13. Total: **143,561 allocated CPU-seconds**.
Local export wall time is recorded separately; no measured CPU time is inferred
from it. The final admission receipt pins the result, reviews and rendered files.
