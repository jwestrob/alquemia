# Tighter local-correlation controls: admitted on the fixed18 panel

The independent saved-only review admits **18/18 native roles (15 new +3 reused) and all35 unchanged checks**. It performed 117 raw/provenance/arithmetic checks. No molecular calculation was run. [Actual review](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_reference_uncertainty_v1/tighter_controls_v1/independent_actual_review_v2/REVIEW.json) and [full evidence](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_reference_uncertainty_v1/tighter_controls_v1/independent_actual_review_v2/EVIDENCE.json) preserve every raw token, method, basis/ghost/electron, residual and cost authority.

## Accuracy of this local approximation on the consumed controls

Errors below are tighter DLPNO-CCSD(T1) minus the original canonical reference, kcal/mol; paired means Dy−La. The original .25 limits remain unchanged.

| Donor | La CP error | Dy CP error | Paired CP error |
|---|---:|---:|---:|
| water | -0.003069 | 0.012664 | 0.015734 |
| acetate_bi | 0.040803 | 0.108710 | 0.067907 |
| acetamide | 0.070324 | 0.160365 | 0.090041 |

Donor-minus-water paired-exchange errors are acetate **0.052174** and acetamide **0.074307** kcal/mol. All18 RHF differences are exactly zero at the admitted printed precision, and all six parent AuxC checks pass. The previous Tight recipe passed 32/35; its failed checks and raw result remain preserved. This supports using the specified tighter local approximation for a separately declared extension of this control comparison. It does not establish physical accuracy, qualify W8/W9 collective energies, or authorize a cluster release.

## What was retained

The recipe is TCutPairs1e−6, TCutPNO1e−8, TCutMKN1e−4 and TCutDO.005, with iterativeT1/CC stopping criteria retained. No change was made to geometry, ghost masks, restricted effective states, ECP46/55, original orbital/AuxC assets or NoFrozenCore. Actual triples-specific thresholds remain separately recorded; this is not a claim that every triples setting was tightened. Maximum terminal CC residual is 9e-09, maximum absolute terminal triples energy change 5.32e-09 Eh; both satisfy their frozen admissions.

Native Energy Check stops remain explicit. **18/18 rows retain one or more failed printed SCF subcriteria**, under the original accepted native stopping convention. These are not recast as complete density/DIIS convergence or state-stability proof. T1 diagnostics span 0.001366902–0.022393821; the original localized-orbital/no-recanonicalization warnings remain in evidence. No new diagnostic threshold or state claim is introduced.

## Actual cost and preserved failures

Only parent Slurm allocations are counted; batch/MPI step rows are retained but not added again. New15 molecular jobs cost **181,224 allocated CPU-s**. Reused three molecular jobs cost **51,800**, reported separately. The cancelled predecessor18518 cost **59,512**; it is not replaced or hidden. The four tighter collectors cost14 CPU-s (including cancelled18807 at0 and failed18519 at2). Total tighter-branch allocation is **292,550 CPU-s**. Original Tight18 plus its collector cost **64,362 CPU-s**, retained as a separate earlier experiment. Canonical reference-production costs are shared prior authority outside this accounting subtree. [Frozen accounting](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_reference_uncertainty_v1/tighter_controls_v1/independent_actual_costs_v1/ACCOUNTING.json) preserves all42 parent rows and the command/output.

The initial independent reviewer used the full-precision FINAL token in a strict comparison to the producer's nine-decimal E(CCSD(T)) token. It stopped on rounding before issuing a verdict. The preserved v1 failure and versioned v2 source correct only that token selection: independent Decimal cycles now use the exact admitted E(CCSD(T)) token, while FINAL/CCSD/triples closure checks remain. No scientific threshold or producer value changed. Successful saved review cost 14.52s wall/13.76 CPU-s on one login thread; the failed reviewer did not record its elapsed cost.
