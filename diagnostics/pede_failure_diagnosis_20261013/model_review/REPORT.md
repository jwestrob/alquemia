# PedE: the current evidence points to structural transfer plus an unqualified response model

**Update after completed saved-component job3341:** the leading diagnostic is now the **solvent response at the saved origins**, before additional DFT on accommodated candidates. The 21.106215 kcal/mol origin spread between two-Ca samples3 and1 consists of −0.150413 from OMOL and +21.256628 from the solvent correction. See the addendum below; this supersedes the next-test ordering at the end of the original review.

Read-only modeling review, 2026-10-13. **The two wrong PedE calls are not explained by a demonstrated execution defect. Their origins already lie between the reference classes, and accommodation shifts them into the La region.** The same La-favoring response also occurs in correct PedE folds. We therefore cannot assign the failures uniquely to folding or to the energy model. The immediate missing evidence is a matched electronic-response check on these actual PedE geometries, alongside the separate saved-structure audit.

This review made no molecular calls, scientific imports, tests or Slurm submissions. It read existing outputs and the exact executed code. [EVIDENCE.json](EVIDENCE.json) pins the files, selected worker modules and calibration source identities. Historical results below retain their original protocols.

## What is actually being scored

The current expression is

`E_M(q) = E_OMOL,vac,M(q) + E_GFN2,ALPB-water,M(q) − E_GFN2,vac,M(q)`.

MACE energies in eV and GFN2 energies in hartree are converted once to kcal/mol. `R = E_Ca − E_La`; larger values are more La-like on this calibrated model scale. These raw totals are neither affinities nor binding free energies. The current GFN2 implementation is **standalone xTB 6.7.1**, with **tblite 0.7.0** checking each ALPB-minus-vacuum component under the declared component-v2 gates. It is not the historical native-ORCA implementation. The original absolute-energy parity gate remains a separate failed historical test.

The executed CPU worker obtains the **total vacuum OMOL energy of every atom in the 192-atom capped context**, using the pinned OMOL extra-large checkpoint, float64, graph charge and physical multiplicity. This is not the older masked atomic-energy descriptor. There is no added dispersion correction, atomic-energy mask, full-protein MM energy or outer-scaffold penalty. Included fragments contribute their local strain and interactions through MACE; omitted protein atoms do not contribute an energy merely because they exist in the physical mapping.

The context is the frozen canonical La-source 4.3 Å fragment union. It contains the same 14 caps and chemical-state signature across the five new PedE sources. Ca/La charges are −3/−2 with multiplicity 1; the dry PQQ state is C14H3N2O8, formal charge −3. Source and synthetic waters are excluded. Both native Ca ions remain in the complete source; only the PQQ-associated ion is exchanged in the local endpoints. The distant spectator is explicitly outside the context and contributes no direct global electrostatic or mechanical energy. Its structural effect can enter through the fold. These state choices are fixed assumptions, not newly demonstrated causes of the two failures. See [preparation results](../../ca_coverage_recovery_20261012/spectator_preparation/FINAL_REPORT.md).

## What can move, and what cannot

At the matched origin, the selector alternates differential Ca/La loads with individual loads, normalizes by physical heavy-atom displacement, and removes redundant directions. It chooses four common angular modes for both metals. The eligible kinds are side-chain torsions and mapped peptide crankshafts. Each metal separately minimizes its **vacuum OMOL proposal energy**, with SLSQP, one origin start, at most 200 iterations, ±0.8 radian bounds and a 0.8 Å physical-heavy displacement guard.

Unselected physical atoms remain fixed. No metal translation, PQQ pose, water or collective-scaffold mode is added by this search. Caps follow source-derived mappings/Jacobians; they are not independent mechanical coordinates. These checks protect geometry, but do not supply the omitted scaffold's deformation energy.

Both metals then evaluate the same `{origin, Ca proposal, La proposal}` pool with the full composite expression. Operational selection retains the origin for an improvement smaller than 0.1 kcal/mol; mathematical row minima are also saved. All required cells must exist. The selection is therefore fair between metals, but its finite candidates need not be stationary on the solvent-containing surface. A successful optimizer is success on the proposal surface, not proof of a solvent-equilibrated pocket.

## What separates the wrong and correct PedE sources

The [completed component/selection audit](../../ca_coverage_recovery_20261012/two_ca_scoring/PED_E_RECOVERY_FINDINGS.md) provides the actual numbers and raw pins:

| Sample | Origin call | Accommodated call | Ca accommodation work | La accommodation work |
|---|---|---|---:|---:|
| 0 | Ca | Ca | −2.123 | −11.792 |
| 1 | Inconclusive | La, wrong | −3.418 | −12.837 |
| 2 | Inconclusive | La, wrong | −2.102 | −12.305 |
| 3 | Ca | Ca | −1.080 | −11.663 |
| 4 | Ca | Ca | −1.292 | −7.264 |

Works are model kcal/mol, rounded here only; `R_after = R_origin + W_Ca − W_La`. Similar La stabilization occurs in correct samples0/3 and wrong1/2. Sample1's selected candidates are interior. Sample2 selects the interior La proposal for both metals, despite its unselected Ca search reaching the displacement bound; correct sample3 does the same. Neither a boundary stop nor the magnitude of La relaxation uniquely identifies the errors.

All paired-coordinate/state/map checks and all 30 metal/candidate solvent-component checks pass. This excludes the specific failures those checks test; it does not establish electronic-response accuracy. The wrong samples select motions involving Glu221/Asp358, with Asn308 also selected in sample2. Significant omitted **vacuum** force loads remain, but such loads also occur in correct controls. More freedom is a hypothesis to test, not a demonstrated repair.

## Calibration and source conditioning are a real transfer boundary

All 25 designated canonical calibration sources are **La-conditioned folds**, including PedE sample2. The fixed calibration rule takes the largest Ca-class canonical contrast and the smallest La-class contrast, separately for origin and accommodated policies. It was not refitted to these two-Ca sources. The current accommodated gap is 3.472171 kcal/mol; the origin gap is 10.297679 kcal/mol. Different bands belong to different scorers.

The new PedE folds condition on two Ca ions. Their local atom inventory matches the calibration policy, but their source geometry and full-source occupancy hypothesis differ. This is a genuine structural-transfer test, not an implementation inconsistency. The two-Ca interpretation is directly motivated by the ExaA crystal evidence; for PedE it remains a modeled occupancy hypothesis. Geometry/confidence, coordination topology and hydrogen-bond comparisons belong to the parallel saved-structure audit. High fold confidence alone would not validate a metal-specific energy response.

## What the completed DFT work does and does not cover

The [actual native DFT donor experiment](../../accommodation_torsion_20260920/FINAL_DFT_REPORT.md) covers 4MAE and two unlabeled PLM sample0 contexts, using ±0.2 radian extra-Asp torsions. OMOL agreed with all 12 endpoint-work signs and all six differential signs; its largest differential error was 1.420383 kcal/mol. That establishes useful response on those paths. **It did not test PedE, these Glu/Asn/catalytic-Asp combinations, or the present optimized candidates.** The published report explicitly leaves coupled-mode and relaxed-candidate validation open.

Adding the historical native-ORCA GFN2 solvent term worsened the PLM differential agreement; [tighter native stopping](../../strict_donor_response_20260923/REPORT.md) changed contrasts by at most 0.080820 kcal/mol and did not remove the discrepancy. That is evidence against treating composite response as automatically validated by OMOL response. It is not a measurement of today's corrected standalone-xTB solvent component. r2SCAN-3c/CPCM and OMOL plus GFN2/ALPB also differ physically; neither treatment is an oracle.

An earlier [protein-only scaffold-proposal experiment](../../scaffold_restart_round_20260922/REPORT.md) included PedE and harmed broader discrimination. Its omitted direct metal/PQQ coupling and absent force-generating solvent make it an inadequate test of a coherent scaffold penalty. Reopening that unchanged generic relaxation rule is not justified.

## Most informative next test

Finish the parallel saved-structure comparison first. If it finds no preparation or gross fold defect, compare the **existing origin and admitted proposal geometries** with an independent electronic treatment. A minimal declared development test can use wrong sample1 and correct sample3, both metals, and all three current pool candidates: 12 fixed-geometry endpoints, with exact deduplication and compatible reuse only. This contrasts a wrong and correct source with similar La accommodation without introducing more optimization or changing bands. A finite manifest must freeze method/state/solvent settings before execution; this review does not launch it.

Compare each metal's work and `W_Ca−W_La`, not absolute cross-method totals. Disagreement would locate a concrete local energy/solvent-response problem. Agreement would shift attention toward source-state/calibration transfer, omitted scaffold/environment and basin sampling; it would not prove that AF3 is wrong. Even a same-context DFT check cannot validate the absent whole-protein mechanical contribution. Do not tune springs, waters, mode counts or thresholds to restore PedE's expected label.

## Addendum: solvent response is the first target after job3341

The [actual completed saved-energy decomposition](../../../workspaces/pede_failure_diagnosis_20261013/energetics_v1/execution/REPORT.md) now identifies where the **origin** differences enter. From correct two-Ca sample3 to wrong sample1, the total contrast changes +21.106215 kcal/mol: OMOL contributes −0.150413 and GFN2(ALPB−vacuum) contributes +21.256628. Relative to canonical La-conditioned PedE, sample1 is displaced +9.128074 overall, made of −2.513729 OMOL and +11.641803 solvent; sample2 is displaced +4.668748, made of −0.721848 OMOL and +5.390595 solvent. These are actual Slurm-analyzed saved values, not a new evaluation in this review.

This substantially narrows the investigation. The large starting spread is already present before the proposed motions, so a new 12-endpoint accommodated-path DFT comparison is **not the first test to run**. The corrected standalone-xTB and independent tblite calculations agree numerically under component-v2; agreement between two implementations of the same method does not establish that its charge response or continuum model is physically accurate for these capped, charged pockets.

Recommended order:

1. **Read the already-computed origin populations and SCF evidence.** Compare both metals, vacuum and ALPB, across all five new sources and the four exact-method La-conditioned controls. Use actual source/core mappings to group target-metal, PQQ and donor-residue populations; retain caps and other fragments explicitly. Examine per-metal solvent components as well as their difference. This can identify charge redistribution, different converged electronic solutions, or sensitivity associated with geometry/exposure. Population changes are method-dependent descriptors, not independently measured oxidation states or a proof of error. A missing raw component remains missing.
2. **If the saved audit reveals no concrete numerical/state mismatch, qualify solvent response at fixed origins with an independent treatment.** A minimal declared diagnostic is the already-identified sample1/sample3 origin pair, both metals, with matched vacuum and continuum calculations: at most eight new endpoints before exact reuse. Select one defensible method/basis/ECP and coherent solvent convention in advance; compare each metal's transfer term and the double difference across the two origins. This tests the dominant observed discrepancy without introducing new geometry or optimization. These consumed, outcome-inspected sources are development diagnostics. Cross-method disagreement still requires interpretation of cavity and electronic-method differences, not automatic acceptance of DFT as truth.
3. Return to the larger accommodated-candidate comparison only if origin solvent response is credible or a specific residual response question remains. Do not fix the current calls by deleting solvent, adjusting radii, altering charge/protonation, or moving reference bands after inspection.

The missing scaffold and source-conditioning limitations remain real, but the new evidence gives a more direct lead than an immediate mechanics rewrite. Geometry differences can cause the solvent sensitivity; the decomposition does not prove that either the fold or ALPB is wrong. The parallel geometry audit supplies the structural evidence separately.

### Existing charge artifacts; no new calculations needed to inspect them

[ORIGIN_CHARGE_ARTIFACTS.json](ORIGIN_CHARGE_ARTIFACTS.json) inventories all **36 existing origin input cells**: nine supported sources × two metals × two media. It pins the actual raw population file, stdout, execution receipt, XYZ and parsed result. It does not calculate populations or group sums.

- The saved energetic result exposes exact cell pins under `rows[].component_input_pins["Ca__origin"/"La__origin"]["vacuum"/"water"]`.
- Each pinned point `RESULT.json` contains `backends.xtb.charges_e`, `backends.tblite.charges_e`, and each backend's `solver_evidence` and `execution` receipt.
- Primary raw populations are in `charges` beside the xTB execution receipt; independent populations are in `scientific.json` beside the tblite receipt. Use the **actual backend execution path**, because reused control cells can live in a different workspace.
- The solvent manifest's `cases[].maps[metal]` identifies `context.core_links` from local atom order to physical source atoms or caps. `context.source_atom_metadata` is in physical-source order, not the charge-array order. `context_preparation.json` and the XYZ pins provide the remaining chemistry/order audit. Never sum populations using unverified positional correspondence.

Source for the addendum is job3341's immutable `energetics_v1/execution/RESULT.json`, pinned by the charge inventory. No new scientific execution, fit, tests or optimizer calls were performed for this addendum.
