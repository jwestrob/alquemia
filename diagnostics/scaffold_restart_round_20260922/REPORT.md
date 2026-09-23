# Nikasha: native restart and collective-pocket result

**Native continuation removes a demonstrated numerical artifact. The tested collective protein-relaxation proposals do not improve the classifier and should not be promoted.** Production fast-PQQ, explicit DFT, reference bands and historical outputs remain unchanged.

## Main findings

The eight confirmed native restarts remove the Q88La 4.805 kcal/mol discontinuity between nearly identical geometries: self/cross transfer differences fall below0.000090 kcal/mol. The first eight seed-only attempts actually used SAD and are retained as failed restart activations. Matching archived GBW plus xtbw activates native continuation. This supports insufficient settling of the original calculation, not a proven unique internal stopping defect or universal ground-state convergence. OnlyLa was tested; no new La/Ca classification follows yet. See [restart report](../native_xtb_restart_20260922/REPORT.md).

For the scaffold test, all24 corrected proposals and all48MACE/96GFN2 scoring endpoints complete. Both metals access the same expanded geometry pool. The old adaptive-band transfer drops from7correct/1wrong to3correct/5wrong; these are transferred bands, **not** a newly calibrated accuracy estimate. A scale change alone does not explain the limitation: both A0AC La-associated folds now rank below the Q9 and Q88 Ca controls, reversing their prior order. No threshold is fitted on these eight consumed sources.

| Source | ΔR versus adaptive, model kcal/mol | Selected Ca / La geometry | Adaptive-band call → candidate |
|---|---:|---|---|
| 1H4I | -9.921453 | scaffold_La / scaffold_La | Ca-supported → Ca-supported |
| 4MAE | -16.681522 | scaffold_origin / scaffold_origin | La-supported → Ca-supported |
| q9z4j7-pqq-la_model | -21.099945 | scaffold_La / scaffold_La | Ca-supported → Ca-supported |
| q88jh5-pqq-la_model | -17.760085 | scaffold_Ca / scaffold_Ca | Ca-supported → Ca-supported |
| a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1 | -8.004827 | scaffold_Ca / scaffold_La | Ca-supported → Ca-supported |
| a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-3 | -10.545813 | scaffold_origin / scaffold_La | La-supported → Ca-supported |
| a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4 | -42.738685 | scaffold_origin / scaffold_origin | La-supported → Ca-supported |
| a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-4 | -57.717670 | scaffold_origin / scaffold_origin | La-supported → Ca-supported |

The A0A3 pair spread contracts8.20→5.66model-kcal. The A0AC Ca/La-conditioned pair contracts15.04→0.057, but both move toward the wrong relative ordering. Agreement between folds is therefore not sufficient evidence of useful discrimination. The still-wrong A0A3 Ca-sample1 is not repaired under transferred bands. All eight pooled contrasts shift downward, by8.00–57.72model-kcal; both MACE and solvent contributions are retained in the full report.

## What physically changed

The source protein-only ff19SB model generated coupled pocket motions with fixed original/Ca-adaptive/La-adaptive donor targets. Source atoms, H inventory, covalent bond lengths, real stereochemistry, peptide cis/trans, fixed exterior, context membership, metal/PQQ coordinates and cap mapping are retained. All proposed heavy motions stay within0.8Å plus the declared numerical tolerance. ff19SB energy is **not** added to the scoring expression.

Every final search reaches the movement guard with significant residual force; none is declared stationary. Limiting moving atoms are mostly outside the compact scoring context. The free-heavy response directions at adaptive targets have cosines0.97066–0.99893 against the original-target control, showing that generic protein relaxation dominates this proposal rule. The original-target response alone already shifts all eight paired contrasts downward. Target-minus-control responses are separately reported, not treated as an exact correction.

This model omits direct protein–metal/PQQ force coupling and solvent in its proposal objective. Donor constraints and overlap guards do not replace those interactions. The observed boundary stops therefore do not establish that physically coherent scaffold accommodation is exhausted. They do establish that this protein-only, fixed-target proposal rule has not earned routine use. The known nativeGFN2 continuity limitation also remains in this unchanged scoring recipe; the restart diagnostic was not opportunistically applied to selected scaffold cells.

Technical recoveries were recorded before scoring: CUDA failed before any force; OpenCL double executed the same model. An initial generic carbon-volume guard incorrectly treated achiral CH2/CH3 atoms as stereocenters. Its942force calls are preserved but unscored. The corrected version protects actual source non-GlyCA andIle/ThrCB centers and their verified substituents. All24 searches were then repeated uniformly. See [mechanics report](../collective_scaffold_20260922/REPORT.md).

## Implementation, validation and cost

| Component | Current status |
|---|---|
| Released fastPQQ / explicitDFT | Unchanged |
| `native_GFN2_Q88_La_xtbw_self_cross_restart_v1` | Diagnostic only;8confirmed restarts |
| `source_ff19SB_collective_matched_donor_target_proposal_v2` |24feasible,0stationary proposals; no FF score term |
| `nikasha_collective_scaffold_geometry_native_OMOL_GFN2_ALPB_v1` | Working opt-in Cartesian common-pool scorer;8/8complete; not promoted |

The shared adapter adds Cartesian receipt admission without changing historical angular candidates. It verifies full parent/target/actual force receipts, source/cap reconstruction and paired state before exact-coordinate deduplication. Required-cell failures remain unavailable; prior scores are separately named. The comparator retains matched-target and q0-control responses, component works, operational/mathematical minima, old-band transfers and the two prescribed structural pair spreads. [Formulae and assumptions](SCORING_METHOD.md); [runnable operations](COMMANDS.md).

Real-artifact checks:40 existing shared parser/pool checks passed;6 corrected collective-admission/algebra checks,11 mechanics checks,7 parent-preparation checks and7 native-restart checks passed. The first broad unittest discovery also encountered an unrelated PLM-export module requiring absent pytest; that module was unrun, not counted as a pass. Parser/algebra checks are separate from the actual molecular integrations above. Explicitly corrupted copies of real fixtures test rejection; no scientific energies were fabricated.

All seven allocations are terminal. Total **50,016 allocated core-seconds and959 requested GPU-seconds**, including failed restart activation, zero-force CUDA failure and superseded mechanics. Actual work:48MACE,112nativeGFN2 attempts (96scoring plus16restart attempts),3,169parent force calls,zeroDFT. Corrected mechanics alone took509s ononeH200/32CPUs; scoring48MACE took36allocation-seconds ononeH200/32CPUs and96GFN2 took249s on64CPUs. These are batch development measurements with different CPU hardware, not matched production per-site throughput. Local preparation/tests/reporting are additional and not fully timed. Full accounting/memory limitations are in [COST.json](COST.json).

## Decision

**Close this generic scaffold-proposal branch; prioritize the demonstrated electronic-continuation remedy.** All80 required saved seed pairs exist for a proposed fixed-pool follow-up on1H4I,4MAE,Q88JH5 andA0A3 Ca-sample1. That would reuse MACE, preserve geometries and evaluate both metals/media consistently. It remains a proposed next analysis; no further restart score, threshold refit or wider campaign ran. [Seed inventory](../native_xtb_restart_20260922/SEED_AVAILABILITY.md).

Jacob separately requested a LanM within-lanthanide investigation while this round finished. Water_basins owns that follow-up; its scope/results are recorded separately. It does not alter the completed La/Ca experiment or turn historic LanM site vectors into independent affinity labels. The larger discriminator-improvement goal remains open.

Full unrounded paired results: [COMPARISON_v1.json](../../workspaces/collective_scaffold_20260922/COMPARISON_v1.json). Compact pinned outcome: [RESULT.json](RESULT.json). No push or default change occurred.
