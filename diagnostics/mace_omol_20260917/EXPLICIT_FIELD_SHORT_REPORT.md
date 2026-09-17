# Explicit field + local MACE: representation check passes, candidate fails

**Retain the production baseline.** This candidate completes cheaply but does
not earn expansion or calibration. All eight new short-component evaluations
and eight native potential evaluations succeeded. No DFT, charge fit, whole
MACE, GB solve or training was rerun. Source states and scientific gates stayed
fixed. The active research goal remains incomplete.

## Actual result

| Representation | Exact direct R | GB R | Short context R | Alpha minus GGR extended / connected |
|---|---:|---:|---:|---:|
| GGR extended | -0.133921 | -54.631776 | -8.104589 | — |
| GGR connected | 10.198140 | -50.267391 | -2.989014 | — |
| Alpha 1F6S | 132.642940 | -180.639855 | -4.115505 | -2.523559 / +1.615627 |
| Alpha 6IP9 | 124.623220 | -173.931139 | -4.284386 | -7.739381 / -3.600195 |

Values are kcal/mol-scale Ca-minus-La contrasts; the last column also includes
native vacuum DFT. Larger relative contrast was predeclared as more La-like.
Only **1/4** ordering comparisons passes. These are two consumed biological
groups, with GGR direct same-assay evidence and qualified cross-study alpha
evidence. They are not four independent biological observations or blind tests.

GGR connected-minus-extended decomposition:

| Contribution | kcal/mol-scale |
|---|---:|
| Vacuum DFT | -23.951208420 |
| Exact permanent-field coupling | +10.332061537 |
| Full physical-boundary reaction field | +4.364385776 |
| MACE short full-minus-core | +5.115575103 |
| **Total** | **-4.139186004** |

The total fails the unchanged absolute 2 kcal/mol partition tolerance. The
whole short term cancels across GGR partitions, as its physical full input is
identical. The cap/core/classical representation changes remain approximate;
this component audit does not identify a unique cause of the residual.

## Transferable findings

The projected CHELPG charge representation passes a stronger, relevant check:
its differential direct coupling differs from the native density result by
-0.099381, -0.541941, -0.455469 and -0.686763 kcal/mol for GGR extended,
GGR connected, alpha1F6S and alpha6IP9. Its partition error is -0.442560.
All **5/5** frozen 1 kcal/mol checks pass. The nearest probes are 1.167–1.202A
from an actual QM nucleus or cap, so this is not only a distant surface test.
It supports the representation on these states; it does not validate all
reaction-field energies or all proteins. Charge-fit refinement alone is not
supported as the main remedy for the observed failures.

Explicit direct coupling removes the previous order-120 kcal mismatch, but
this comparison changes both the learned context and the coupling expression.
Do not attribute the whole difference uniquely to one learned component.
The remaining few-kcal partition/ordering errors are material. The two alpha
structures differ in water inventory (2 versus 3); their difference cannot be
called pure structural noise or repaired by selecting the favorable structure.

## Model and limitations

Protocol `normalized_vacuum_density_ff19SB_full_OBC2_POLAR_short_hybrid_v1`:

`A_M = E_DFT,vac(core,M) + sum(q_env * phi_vac,M) + G_full(Pq_QM,M + q_env) + T_short(full,M) - T_short(core,M)`.

`R = A_Ca - A_La`. Native quantum Hartree and MACE eV are converted once with
project constants 627.509474 and 23.06054783061903. Direct coupling uses native
atomic-unit potentials at the exact environmental charge positions. The point
comparison uses the recorded project Coulomb constant; its tiny convention
difference is retained in the reported comparison, not fitted away.

Full GB, boundary charges, source states and six whole short components are
reused exactly. There is no OMOL total, MACE electrostatic/electron scalar,
core CPCM, added QM internal Coulomb, relaxation or entropy. The density is
frozen in vacuum: electronic response to the protein is absent. The GB term
uses its projected monopoles, not an exact density functional. The learned
short/field decomposition is not a unique physical separation.

Absolute reference, calibrated classification, combined gradient and response
correction remain null. New force files differentiate only the local learned
component. This is not a binding-free-energy result.

## Verification and cost

Manifest SHA `90d553d2480a1efb3911fdc9169774610c9b124d40db7e17e6b4ad45af32d439`.
Four source/preparation/cache/partial-failure tests passed in10.525s; the
actual-result test was initially skipped, then passed in0.031s after execution.
Three existing short-engine receipt/component tests passed in1.138s. Eight
distinct tests pass; no unavailable executable was replaced by fabricated
output. Both preparation and frozen MACE-environment dry-run passed.
All4 component algebra checks pass (maximum1.17e-10kcal); inherited full-GB and
short-engine numerical qualifications remain explicitly separate.

- Potential job1200980:29s ×8CPUs =232allocatedcore-s;93.964actualCPU-s.
- Short job1200981:63s ×16CPUs =1008allocatedcore-s;83.743actualCPU-s;
  **63GPU-allocation-s**. Eight evaluations sum7.561124model-s.
- Total incremental allocations: **1240core-s**,177.707actualCPU-s;zero failures.
- Native utility time sums88.160987s across eight concurrent processes;
  peak individual RSS316196KiB. Short peak RSS1640488KiB;
  peak allocated GPU480389120bytes, reserved566231040bytes.

Local prepare/test/collect resource receipts are retained separately. Archived
DFT, charge fitting, whole inference and GB costs are additional; this is not
an end-to-end production cost measurement or matched speedup claim.

**Three separate judgments:** numerical execution and direct representation
checks pass; predictive/partition usefulness fails; incremental execution is
inexpensive, while a full production cost remains unmeasured for this model.
Recommendation: retain baseline and reject this frozen candidate for expansion.

Full results/cost/source receipts are under
`workspaces/mace_omol_20260917/explicit_field_short_v1/`,
`explicit_field_short_result_v1.json`, and `explicit_field_short_cost_v1.json`.
Compact unrounded result: [EXPLICIT_FIELD_SHORT_RESULT.json](EXPLICIT_FIELD_SHORT_RESULT.json).
Read-only next command and exact replay operations:
[EXPLICIT_FIELD_SHORT_COMMANDS.md](EXPLICIT_FIELD_SHORT_COMMANDS.md).
