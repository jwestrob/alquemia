# Water-reference decision tree

**The Dy-only fit stays held.** This decision tree interprets the declared48-cell water test; it adds no calculations, parameters, tolerances or validation observations. Its purpose is to choose the next physical-model decision, not extend an audit indefinitely.

Sources: [PLAN](../../PLAN.md), [CURRENT](../../CURRENT.md), [water48 specification](../../reference_chemistry_v1/water_force_qualification_v1/SPEC.md), [La baseline audit](../../reference_chemistry_v1/la_baseline_audit_v1/REPORT.md), [saved training controls](TRAINING_REPORT.md), and [La force localization](LA_DIAGNOSIS_REPORT.md). Reserved24 outcomes remain unopened.

## 1. Establish what the calculation actually resolves

Require the complete declared native denominator and unchanged state, ghost/ECP, basis, SCF/CC and termination admissions. An unavailable or failed cell leaves its associated derivative unavailable; it is not evidence against either force field. Reuse completed cells. Any recovery must address a demonstrated cause without silently changing the scientific recipe.

For each individual metal **and separately Dy−La**, use the frozen0.1kcal/mol/Å two-step criterion. Summarize h=0.005Å only when the0.005/0.010Å estimates agree at that scale. Otherwise retain both and mark the quantity unresolved: no third step, extrapolation or fitted tolerance follows automatically. This resolution is not a physical-accuracy error bar.

The coordinate derivative is −sum(water forces)·u, where positive displacement moves the intact water away from the fixed metal. Compare the saved PBE0-D4 analytic derivative with the **uncorrected PBE0-D4 complex** finite difference. CP derivatives are different observables and must not be substituted for that numerical check. Preserve D4 separately and central EnGrad/displaced-SP final-integration/COSX signatures.

## 2. Branch on the actual evidence

| Observation | Interpretation | Next action; what remains blocked |
| --- | --- | --- |
| Two-step estimates are unresolved | Truncation, numerical noise or state/energy-stage changes remain possible. | Keep that derivative unavailable. Inspect the already-preserved native stages to identify a concrete cause; do not fit, discard the inconvenient step or infer model error. |
| Resolved PBE0-D4 complex finite differences disagree with the saved analytic derivative | The sampled energy/force interface is inconsistent at the declared scale. This can involve EnGrad/SP stage differences, not necessarily a bad analytic-gradient formula. | Keep those analytic targets ineligible. Correct only an identified representation/interface cause through a separately declared version; no force-field repair or energy/force mixing to hide it. Other admitted method energies remain valid observations. |
| Analytic consistency passes, but CP substantially changes the radial response | Basis sharing affects the target; the uncorrected force is not a CP force. | Keep CP and uncorrected surfaces separate. Use existing matched CC/PBE0 and donor-CP evidence to choose a common target before fitting. Never combine CP energies with uncorrected gradients. |
| AMOEBA La agrees with resolved CC-CP response but misses uncorrected PBE0 | The original La failure may principally reflect the chosen reference target on this motion. | Retain published La unchanged while assessing Dy against the same target and inspecting the already-declared donor comparison. Do not fit La toward PBE0 merely to clear its old control. This one water result does not clear the donor or protein gates. |
| Resolved methods disagree materially, with no supported target choice | Method sensitivity remains consequential. Agreement with whichever result is favorable is not validation. | Hold fitting and record the target uncertainty. Use the existing correlated donor and electronic-state evidence to make one explicit target decision; no method selection by metal, observable or desired preference. |
| AMOEBA misses a radial response consistently across the credible, resolved references | A local radial model discrepancy is supported on this fixed fixture. Attribution to a particular interaction term is still unproven. | If only Dy fails and La passes the chosen target, reconsider only the declared Dy-water stage after reference eligibility review. If La also fails, the unchanged fixed-La/Dy-only family cannot qualify the pair; make a separate paired-model repair decision before fitting. No additional parameters are specified here. |
| Both metals agree with the selected reference on this motion | Water radial repair is not supported by this test. | Keep water parameters fixed; proceed to the already-prepared donor/transfer qualification. Do not repeat water fitting or claim hydration ensembles, full-vector forces or protein selectivity are validated. |

Multiple rows can apply. Resolve numerical eligibility first, then distinguish reference-target sensitivity from model discrepancy. Small Dy−La error with large individual errors is cancellation, not clearance of the La or Dy control; it can still distort either solvent ensemble.

## 3. Preserve a coherent pair and take one next step

The published La model used MP2 geometries, CP-DLPNO-CCSD(T)/ZORA interaction curves and hydration evidence. Water48 uses canonical CCSD(T) or PBE0 with the shared frozen-f ECP/basis convention. Its CC result is therefore **not an exact reproduction of the published calibration or physical truth**. Basis, relativistic/core, Dy-state and ensemble uncertainties remain. Native stopping alone does not resolve them.

Choose one compatible energy/force target for both metals, donor environments and thermodynamic legs. Water48 supplies only a sampled radial derivative: it cannot replace full CP force vectors in the existing rigid-body training objective. If a CP target is chosen, compatible force targets or a separately declared objective are required before any fit. First use the ongoing donor-CP comparison to decide whether that additional work is warranted.

At the outcome review, make one recorded decision: repair a demonstrated numerical interface, retain La and qualify Dy on the selected common target, or reject the fixed-La assumption and declare a paired repair. Preserve all failed controls. None releases production sampling: candidate numerical forces, donor/protein transfer, unopened validation, periodic polarization and ensemble/free-energy gates in PLAN still apply. Existing topology and alchemical-path correctness remain useful without being promoted to physical qualification.
