# Prospective paired-donor development

**Proposal only: no fitter, new energies, candidate selection or validation opening.** The held seven-parameter Dy-only protocol remains unchanged. First complete the separately frozen water-candidate derivative/transfer test. Job17576 selected water parameters `[3.710148544089586, 0.2001607792689767, 3.389554168967775, 0.16712500000000005]` and passed its four consumed-center controls; this does not establish transfer to hydrated donors or protein.

## Conditional question and minimal family

If the water candidate passes its numerical and consumed water-radial controls but leaves donor-specific errors, ask whether **eight pair parameters** can represent both metals' donor responses on the same uncorrected PBE0-D4 surface. Freeze all four selected water parameters and their epsilons. Allow only radius and damping for each metal's carboxylate oxygen and amide oxygen. Mono- and bidentate acetate share one carboxylate pair; no environment-specific values.

| Pair | Initial Rmin, Å | Initial damping | Fixed epsilon, kcal/mol | Small-system VDWPR class / POLPAIR type |
|---|---:|---:|---:|---|
| La–carboxylate O | 3.6318 | 0.1518 | 0.2280 | 403–500 /403–500 |
| La–amide O | 3.3000 | 0.2987 | 0.3019 | 55–500 /105–500 |
| Dy–carboxylate O | 3.347237658604727 | 0.139 | 0.3795688850865599 | 403–501 /403–501 |
| Dy–amide O | 3.2657502966294456 | 0.139 | 0.3795688850865599 | 55–501 /105–501 |

Use the prior declared radius bounds, 0.85–1.15 times each **resolved pair** initial Rmin, and damping bounds0.05–0.50. These are exploration limits, not uncertainty intervals. La values come from explicit SI overrides; Dy values come from native CUBIC-MEAN/HHG resolution, not atomic diameter. Preserve the La amide-C override `(class53,500): Rmin3.9541, epsilon0.8295`; preserve Dy amide-C native mixing. Freeze all other ion/ligand parameters, epsilons, permanent multipoles, intrinsic polarizabilities, bonded terms, water model, geometry and polarization settings. The independent audit found no baseline implementation defect; changing La is a new reference-specific hypothesis, not a correction to the published calibration.

## Exact consumed rows

Use the existing complete28-row dataset, with no new geometry or reference calls:

| IDs (M=La,Dy; F=acetate_mono,acetate_bi,acetamide) | Rows | Role |
|---|---:|---|
| `water_central_M`, `aquo8_central_M` | 4 | Reuse frozen candidate controls and matched aquo energy anchor |
| `hydrated_F_{minus,central,plus}_M` | 18 | All hydrated donor translations at−0.15/0/+0.15 Å |
| `matched_dry_F_central_M` | 6 | Exact water-removed controls, unchanged surviving coordinates |

This is **three donor chemistries**, six wet/dry categories plus two water categories; it is not eight independent donor chemistries. The24 donor rows form12 matched pairs. Other previously inspected dimer, aquo9 and CC/CP evidence remain contextual diagnostics; they are not silently added as fitting targets. The prospective transfer run supplies the new-water baseline on all18 hydrated rows; six dry rows are unchanged. Do not infer those post-water donor errors from the original baseline.

Dataset: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/matched_dry_v1/training_v1/RESULT.json`, SHA256 `99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb`. It retains exact native/reference/coordinate/type/graph pins for every row. Water selection: `.../paired_water_repair_v1/joint_refinement_v1/run_v1/execution/RESULT.json`, SHA256 `fab7016c4335a82dea983368abc6f23201dd83df0324f7794e5e0c978fd8e190`.

## Symmetric development targets

Retain the existing net-force and geometric-centroid torque/1Å definition for complete covalent bodies. Weight ion/ligand/water categories equally, waters equally within category, and translation/rotation equally for each non-ion body. Preserve full atomic/internal-force residuals as diagnostics; do not fit them into metal parameters.

Proposed loss has four equal blocks: balanced energy, individual La rigid response, individual Dy rigid response, and matched Dy−La rigid response. Force scale remains2kcal/mol/Å and energy scale1kcal/mol. Within the energy block retain the old relative proportions:0.50 aquo-referenced exchange errors across12 donor pairs;0.30 wet-minus-matched-dry double differences across3 families;0.20 radial-work errors, **equally including both metals** across12 endpoint works. Use each block's mean square. These correlated cycles and force components are deterministic engineering weights, not independent observations or a statistical likelihood.

For D(q)=E_Dy(q)−E_La(q), fit only residuals of C(q)=D(q)−D(aquo8), water-addition D(wet)−D(dry), and same-metal E(±)−E(0). No raw unlike-element or unlike-composition energy offsets, absolute donor binding claim or affinity target. Reuse the frozen candidate aquo anchor. Report paired radial works separately as well.

Retain existing transfer resolutions: every energy error≤1kcal/mol in magnitude; individual La, Dy and paired rigid RMS≤2 and maximum body component≤5kcal/mol/Å, with every fixture/family retained. No scalar-loss improvement can hide an individual failed control. Water controls must remain unchanged. Freeze an executable objective, finite joint optimizer budget and exact stopping rules before any fit; this document does not choose a budget or run.

## Identifiability, reference limits and protein mapping

Radius and damping can compensate; a central energy alone cannot identify them. Wet/dry and three radial positions constrain environmental and radial response, but mono-acetate torque/orientation is a harder test of a shared oxygen-pair model. Net forces and torques obey global conservation, so a large residual-vector length does not establish eight independent constraints. Before interpreting a selected candidate, report normalized weighted-Jacobian column norms, rank at relative1e−6, condition (>1e6 flagged) and proximity within1% of bounds. Joint bounded trials address coupled directions; a budget or local stop alone cannot establish family impossibility. No automatic epsilon addition, bound widening or re-fit after validation.

Targets remain uncorrected ECP-PBE0-D4 energies **and gradients**: La46/Dy55 frozen-f cores, effective singlets with physical labelsLa1/Dy6. Native completion is not physical accuracy. The existing coarse acetamide-La0.15Å analytic/secant discrepancy remains flagged (1.713 versus1.485kcal/mol/Å screen); it is not proof of derivative failure and cannot be dropped or loosened. The existing dataset explicitly requires a root reference-eligibility decision before fitting. Water48 small-step qualification does not automatically qualify donor gradients. Published La used CP-DLPNO-CCSD(T)/ZORA and hydration calibration; the consumed canonical donor CP comparison shows donor-dependent reference sensitivity (amide1.401006 and biacetate0.254895kcal/mol PBE0-D4-minus-CC gap magnitudes). Keep these separate; no CP-energy/unCP-force mixture and no claim that a PBE0 fit recovers published physical calibration.

Protein transfer remains a separate hypothesis. `periodic_mapping_v1` record SHA256 `edaa4cc6cbeb8b82670e0314423432667102640eb4cf0df84e2eb0bc300f34e7` maps carboxylate VDWPR to class31 and POLPAIR to types143/161; amide O to class5/types5,11,153; fixed amide-C correction to class3. Fragment-only ACE226/227 recipients must also be covered in their own namespace. Small acetate403 conflicts with BIO water-H403: direct numeric copying is forbidden. Native protein oxygen/carbon multipoles and polarizabilities remain unchanged and differ materially from SI small donors. BIO water Urey–Bradley differs too. A successful small-system fit cannot promote this map as transferable.

## Next decision

- Numerical water-candidate failure: stop this candidate; no donor fitting around it.
- Water radial failure: resolve water representability before permitting donor parameters to compensate.
- Water checks pass, donor checks fail: the eight-parameter hypothesis is a concrete next **development** candidate, subject to the existing reference-eligibility decision and a new frozen execution specification.
- All consumed donor checks already pass: no unnecessary donor fit; proceed only to the separately declared validation decision.

The reserved24 orientation/protein outcomes remain unopened throughout this proposal and any later development fit. No new scoring, protein simulation or physical qualification follows from this document.
