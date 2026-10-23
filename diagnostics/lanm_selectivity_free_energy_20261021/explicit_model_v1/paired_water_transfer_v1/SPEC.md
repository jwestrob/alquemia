# Frozen water candidate: derivative and consumed-configuration transfer test

## Decision and execution boundary

**Question:** Does the single frozen four-parameter water candidate have consistent forces and extend beyond its four fitted centers to already-consumed water radial and hydrated donor configurations?

This specification is written before inspecting pilot 17567 outcomes. It changes neither that pilot nor its selected candidate. Root may release the following finite test only after admitting the actual pilot result, its selected-vector/source/input pins, and complete candidate records. Record the pilot's development, numerical-stop, rank and bound flags separately. A failed/incomplete pilot supplies no executable candidate; a complete but development-failing or budget-limited candidate requires an explicit root diagnostic decision and must retain those labels. No automatic execution follows this file.

Use the candidate's exact La/Dy water-O pair radii and damping, its fixed epsilons, and the same pinned CPU executables, native admission and small-molecule SI parameter file. Preserve charges, covalent terms, permanent multipoles, ion/water polarizabilities, water-H and donor pairs, mutual polarization settings, coordinates, bonds and atom ordering. No parameter selection, fitting, optimization, new reference calculation or reserved validation is included. The protein BIO water model and periodic transfer remain separate.

## Exact inventory and finite cost

Let `M={La,Dy}`, `F={acetate_mono,acetate_bi,acetamide}`, and `S={minus,central,plus}`. These sets define exact IDs, not permission to discover additional cases.

| Group | Exact IDs | Treatment | Maximum new native programs |
|---|---|---|---:|
| Fitted centers | `{water,aquo8}_central_M` (4) | Reuse candidate E/analytic forces; one full Cartesian `testgrad system.xyz Y Y 0.00001` per center | 4 |
| Consumed water radial | `water_{minus,plus}_M` (4) | One `analyze E` and one analytic `testgrad Y N` per cell | 8 |
| Consumed hydrated donor panel | `hydrated_F_S_M` (18) | One analyze and one testgrad each; use full Cartesian finite differences for the six `central` cells, analytic-only for the other twelve | 36 |
| Matched dry controls | `matched_dry_F_central_M` (6) | Reuse unchanged saved E/forces after proving no water-O type/class recipient exists | 0 |

**Ceiling: 48 new native program invocations.** There are ten full Cartesian checks, containing 4/25 atoms at the fitted water/aquo8 centers and 29/26/31 atoms in the respective hydrated monoacetate/biacetate/acetamide centers, for each metal. The native routine evaluates two displaced energies per Cartesian component: **1,380 internal displaced `energy()` evaluations** across the ten checks, in addition to up to 26 analytic gradient evaluations and 22 analyze evaluations. Program counts are not energy-evaluation counts. Exact eligible cache reuse can reduce these ceilings; never repeat an identical completed input.

In the pinned `testgrad.f`, argument `0.00001` is the **total separation**: coordinates are evaluated at ±0.000005 Å and divided by 0.00001 Å. Require all atoms active and all analytic/numerical components present. Keep this convention and the original 0.001 kcal/mol/Å component tolerance unchanged.

The 28-row training result supplies the four fitted centers, eighteen hydrated rows and six dry controls. The four water radial rows come from the earlier consumed panel and its existing ECP-PBE0 results. All original AMOEBA baselines and reference E/forces are reused. Unavailable evidence remains unavailable; it is never replaced by zero or by a new reference call.

## Stage A: numerical force consistency

Check the four fitted centers first. For every full Cartesian check, require finite native energy/gradients, successful polarization/native admission, and maximum absolute analytic-minus-numerical gradient component ≤**0.001 kcal/mol/Å**. At fitted centers, additionally replay the frozen candidate energy within **0.0001 kcal/mol** and analytic force components within **0.001 kcal/mol/Å**, the existing engine comparison tolerances. Compare gradients to gradients; force is minus gradient.

Only after all four pass may the consumed transfer cells run. Apply the same Cartesian gate to all six hydrated central cells. A numerical or native failure stops new scheduling and closes this candidate's transfer qualification. Preserve already-running siblings and all attempted receipts. No changed finite-difference step, tolerance, solver, extra search or retry is included.

## Stage B: frozen transfer comparisons

Admit exact source coordinates and state: charges +3 for water/aquo8/acetamide, +2 for acetate complexes; La physical label 1 and Dy 6; saved ECP references retain cores 46/55 and effective singlets. Compare to the same **uncorrected PBE0-D4** surface used by the pilot. Do not mix CP/CC energies into this force comparison; retain the pilot's separate Water48 method-sensitivity report.

**Path signs must follow the actual saved motion.** In the original `water_minus/plus` fixtures the metal alone moves by q times the stored water-to-metal unit vector, with q=±0.15 Å; positive q increases the metal–water separation. Its derivative is `−F_metal·u`. In hydrated `minus/plus` fixtures the intact donor ligand moves by q times its stored outward direction, with waters/metal fixed; its derivative is `−sum(F_ligand)·u`. Read and verify these motions from the pinned coordinates. Neither path is silently reconstructed as the Water48 rigid-water translation. The 0.30 Å endpoint secant is a coarse shape diagnostic, not the small-step derivative qualification.

Report original baseline and frozen candidate side by side, retaining every individual La/Dy and Dy−La result:

1. **Water and donor radial work:** for each metal and ±0.15 Å state, `E(side)−E(central)`, plus its matched Dy−La difference and reference error. This tests a finite curve beyond the fitted center.
2. **Balanced exchange:** for every water side or hydrated state q, `C(q)=[E_Dy(q)−E_La(q)]−[E_Dy(aquo8)−E_La(aquo8)]`. Each model uses its own matching aquo8 denominator. Both sides of this exchange contain one ion of each element and the same total ligands/waters; no raw unlike-element or free-fragment offset is interpreted.
3. **Matched water addition:** for each central donor family, compare `[E_Dy(wet)−E_La(wet)]−[E_Dy(matched-dry)−E_La(matched-dry)]`. The unchanged dry controls require no calls. This measures a metal-dependent many-body water-addition effect, not a unique metal–water pair contribution.
4. **Rigid response:** use the existing geometric-centroid force/torque helper. At each cell, weight ion, complete donor ligand where present, and water categories equally; average waters within their category, with translation/rotation equally weighted per non-ion body. Report La, Dy and paired category-balanced RMS and maximum body-component errors separately. Retain internal atomic-force residuals as diagnostics.

**Frozen development resolutions:** every listed radial-work, paired radial-work, balanced-exchange and matched-water-addition error must be ≤**1 kcal/mol** in magnitude; every individual and paired per-fixture rigid RMS must be ≤**2 kcal/mol/Å**, with maximum body-component error ≤**5 kcal/mol/Å**. These are the existing development scales applied to explicit transfer checks, not new physical error bars. Do not average failures away across families, metals or radial positions. Report error changes from baseline, including any deterioration, regardless of threshold status. No threshold or weighting changes follow inspection of the candidate.

## Interpretation and resources

| Outcome | Allowed conclusion |
|---|---|
| Candidate native/Cartesian gate fails | Numerical derivative or admission failure; no transfer claim |
| Fitted centers improve, water radial transfer fails | Fit is local to its fitted centers |
| Water radial passes, hydrated donor transfer fails | Limited water-only transfer; coupled donor/water behavior remains inadequate |
| Errors improve but some frozen development checks fail | Partial improvement with explicit remaining failures; no promotion |
| All numerical and consumed-transfer checks pass | Useful **development transfer** on already-consumed configurations; not independent validation, bulk hydration, protein transfer or affinity prediction |

Root owns a separately released allocation, at most eight CPU workers with one native process and one numerical thread each. No GPU or login molecular work is proposed. Record native program receipts, the internal finite-difference count, aggregate child CPU once across the run, parent CPU, elapsed and allocated cost. Reuse unchanged baseline and candidate inputs by exact physical-input/source identity. Preserve a complete cell denominator and stop on the first observed native/nonfinite/polarization failure. No reserved24 outcome is opened, and no execution or parameter change is authorized by this specification alone.

## Pinned existing evidence

All paths below are relative to the project root. Candidate result/key/cache pins remain unassigned until root admits the actual pilot result; no outcome has been read to prepare this specification.

- Training: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/matched_dry_v1/training_v1/RESULT.json`, SHA256 `99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb`.
- Original consumed AMOEBA panel, including water radial inputs: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/panel_v1/RESULT.json`, SHA256 `b256fdca923c0f476651e73250b93672eb49f7634fbceca66d2e6ea1d398a052`.
- Water radial QM result base: `workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/ecp_remaining_v1/tasks/ecp_pbe0__CELL/execution/RESULT.json`:
  - `water_minus_La`: `e889750771ca1db81a72bc02f492df9f82c2968e23bbe70745e7d63cfa2e1206`
  - `water_plus_La`: `404455607dcdfea0f1857f4da68cd6ba3ab4affcf4057275e6bde80c5e3c2ab8`
  - `water_minus_Dy`: `098d03cc2893463fd3c19a340f4bef6dadd0414cc9371d76c80a7e4912a659d9`
  - `water_plus_Dy`: `05dd6ce6305828bb17eae27e6c3712ce8427ae1b13f0c2db459766c6ae58b4f1`.
- Rigid-response helper: `diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/matched_dry_v1/OBJECTIVE.py`, SHA256 `bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819`.
- Native finite-difference convention: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/tinker_source/source/testgrad.f`, SHA256 `5b3523b67e95441bbc0d45582b0055c263cbda13fd0a8b03dc81d8ce33924333`.
- Pilot specification: `diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_repair_v1/SPEC.md`, SHA256 `985977202ec5f1da9b93942c6a2e0d44e7a07a2edac674bd72f7ef76edbc22f3`.
