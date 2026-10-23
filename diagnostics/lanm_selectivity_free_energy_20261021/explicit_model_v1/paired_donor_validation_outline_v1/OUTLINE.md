# Conditional validation after donor development

This is a source-only execution outline. No candidate is bound, reserved energy or force outcome opened, or computation released. Preserve the existing 24-cell reserved manifest unchanged. A passing donor-development fit would justify testing this panel; it would not establish a qualified protein Hamiltonian.

## Frozen inventory and order

| Stage | Configurations, each with matched La/Dy | Cells | Question |
| --- | --- | ---: | --- |
| A | Aquo8 water-H orientations −15°/+15°; hydrated mono-acetate, bi-acetate and acetamide ligand orientations −15°/+15° | 16 | Do fixed water/donor parameters transfer to untrained orientations within the original small-molecule Hamiltonian? |
| B | Exact V3 ACE–Asp3–NME, ACE–Glu11–NME, ACE–Asn5–NME and ACE–Trp9–NME | 8 | Do those pair corrections transfer to protein multipoles, polarizabilities and capped chemical environments? |

All 24 slots are currently geometry-eligible. Stage A contains eight paired configurations; B contains four. Retain every slot even if execution or reference admission fails; missing results are unavailable, never zero or silently excluded. Stage B uses V3 only; withdrawn V2 cannot substitute.

The native reference manifest declares exactly 24 `ecp_pbe0` cells and no free fragments. La/Dy physical multiplicities remain 1/6, separately from the effective singlet frozen-f ECP46/ECP55 calculations. Acetate/Asp/Glu complexes have charge +2; water/amide/Asn/Trp complexes +3. Coordinates, atom order, covalent graph and charges stay fixed.

## Entry and execution gates

1. Freeze one actual terminal donor candidate, its eight donor values, fixed six water values, complete development result, parameters, source and fit decision. Require every original development control and water invariant to pass; a budget stop is not a passing candidate. If development fails, keep this panel sealed.
2. Admit the frozen candidate through declared numerical derivative and exact-input replay checks on consumed configurations. Keep existing numerical tolerances: energy replay 1e−4 kcal/mol and force-component/Cartesian derivative consistency 1e−3 kcal/mol/Å. These check the implementation, not physical accuracy. No fitting during these checks.
3. Before opening reserved outcomes, bind the immutable geometry/native manifests, candidate and a finite execution/scoring specification. Reuse the established scoring conventions and thresholds below. Account separately for native program invocations and internal Cartesian finite-difference energies. Root owns submission.
4. Execute and admit all 16 Stage A cells under unchanged small-system typing, reference recipe and native solver settings, then score the complete stage. Preserve both individual metals and paired results. A failed or unavailable Stage A stops progression to Stage B; it is a recorded validation outcome, not grounds to alter an orientation, threshold or candidate.
5. Only after Stage A passes, finish the explicit protein mapping and reference-scope checks below; execute and admit the eight Stage B cells, then report them as a separate transfer test. Success of either stage does not release periodic dynamics or an exchange-free-energy claim.

The existing geometric gates remain unchanged: internal-distance preservation 1e−9 Å after 12-decimal serialization, immobile unrotated atoms, heavy inter-body contacts ≥1.8 Å and H-containing contacts ≥1.2 Å. Use the frozen axes/signs/angles, with the original axis fallback rule. Protein V3 retains its cap/source/contact/typing evidence; there is no new geometry search.

## Scoring, fixed before reserved readout

Carry forward the established development limits: each declared balanced-energy error ≤1 kcal/mol in magnitude, and each individual La, individual Dy and paired rigid-response fixture RMS ≤2 and maximum component ≤5 kcal/mol/Å. Use the existing category/body weights, net forces and centroid torques divided by 1 Å. Whole-atom/internal force residuals remain separate diagnostics. Passing an aggregate mean cannot replace the per-fixture limits. Do not reconstruct a fit loss or invent radial-work targets for the orientation-only panel.

For Stage A, report eight aquo-referenced contrasts with the unchanged sign convention

`C(g) = [E_Dy(g) − E_La(g)] − [E_Dy(aquo8 central) − E_La(aquo8 central)]`.

Compare model minus reference. For each of the four orientation families, also retain La and Dy plus-minus work and their paired difference; this cancels composition-specific offsets without new free-fragment calls. Their errors use the existing 1 kcal/mol limit. Apply individual/paired rigid criteria to all eight configurations (24 force fixtures), retaining water and ligand contributions separately.

Stage B has four paired configurations and therefore 12 individual/paired rigid fixtures. Treat each intact capped fragment as one ligand body; retain atomic/cap-local residuals diagnostically. There is no protein radial or orientation pair in this panel. For aquo-referenced protein energies, use a matching BIO-water model aquo8 anchor at the same fixed reference geometry; reuse exact admitted anchors if available, otherwise declare the two additional model evaluations before execution. An SI-water/BIO-protein combination must not be mislabeled a single consistent protein-to-water exchange model. No isolated-fragment energy is needed for these paired cycles.

The reserved preparation originally froze geometry, not a numerical accuracy threshold table. This outline explicitly carries existing development tolerances into prospective validation; the execution specification must preserve these declarations before any reserved outcome is read. No post-read threshold adjustment or parameter refit can preserve the status of these cells as independent validation. Any later use for development consumes the entire affected validation set and requires new independent evidence.

## Protein mapping is a substantive transfer

Use the exact BIO parameter file and retain its permanent multipoles, local frames, intrinsic polarizabilities, bonded terms and nonmetal topology. Map operators by chemical role: VDWPR uses classes, POLPAIR uses types. Metal types/classes 900/901 are the separate proposed protein namespace; do not copy the small-system SI key.

- Carboxylate O: class31; Asp type143 and Glu type161. Shared mono/bi donor correction remains shared. No carboxylate-carbon override is invented.
- Amide O: class5; represented backbone type11, Asn side-chain type153 **and ACE cap type227**. Include any other represented class5 oxygen found in the actual V3 inventory. The full-box mapping omitted ACE because those caps were absent there.
- Preserve the La amide-carbon VDWPR on class3, including backbone C9, Asn C152 and ACE C226; keep Dy carbon mixing unchanged. No new cap parameter fit.
- Water recipients in BIO are O type402/class73 and H403/class74. Small-system acetate type403 must never become a BIO water-H POLPAIR. The reserved protein fragments contain no waters; matching BIO aqueous anchors still require the correct water mapping.

Source values differ materially: small carboxylate O polarizability 1.931 Å³ becomes protein 0.837; small amide O 0.914 becomes protein 0.837/0.834. BIO ACE O227 is 0.837 and C226 is 1.334. Preserve these differences. BIO water also has Urey–Bradley equilibrium distance 1.5537 Å versus SI 1.5337 Å; rigid fixed-water internal terms can cancel in paired contrasts, but flexible-water models are not identical. Resolve every represented metal–class pair, polarization recipient and local frame explicitly, and qualify the mapped native numerical implementation before interpreting Stage B.

## Reference and budget limits

Use the existing uncorrected ECP-PBE0-D4 reference recipe and basis/state pins consistently for energies and forces. The passing small-step central acetamide translation and water tests do not qualify every orientation, Cartesian force component or protein environment. Preserve the earlier coarse-screen failure and CC/CP method sensitivities; do not splice CP energies into an uncorrected force surface. A native or derivative-admission failure yields an unavailable or reference-limited comparison, not a model pass or automatic solver retry.

The frozen denominator is 24 reference E/gradient evaluations and 24 candidate model E/force evaluations before exact reuse, staged 16 then 8. Native analyze/testgrad program counts, mapped BIO aquo anchors and numerical derivative checks must be separately declared from that scientific-cell count. No extra free fragments, optimization, MD, fit loops or independent-validation expansion is implied. This outline opens no reserved numerical outputs.
