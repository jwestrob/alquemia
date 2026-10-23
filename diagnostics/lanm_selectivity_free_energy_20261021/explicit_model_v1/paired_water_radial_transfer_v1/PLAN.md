# Frozen radial-water candidate: consumed donor transfer

## Question and release boundary

Does a single frozen candidate from `paired_water_radial_repair_v1` retain numerical force consistency and improve the already-consumed hydrated donor configurations?

This is a source-only prospective plan, written without reading the radial pilot's execution outcome. Candidate result, selected-vector and run-release pins remain unassigned. Root decides whether to execute after reviewing the actual finite radial pilot. No native calls, fitting, new reference data, reserved24 outcomes, dynamics or production changes are authorized by this file.

Admit only a complete, error-free `bounded_radial_refinement_complete` result with its exact source/config, selected eight-cell model, unchanged four parameters/bounds/epsilons and full native provenance. Preserve its stopping reason, selected per-fixture controls, iteration-origin rank diagnostics and bound flags. A complete development-failing or budget-limited candidate requires an explicit root diagnostic decision; a failed or incomplete candidate is not executable. Passing the water fit is not itself a transfer result.

## Inventory, reuse and numerical gates

Retain the original prepared32 geometry/state/reference inventory byte-for-byte. Let M={La,Dy}, F={acetate_mono,acetate_bi,acetamide}, S={minus,central,plus}.

| Group | Exact cells | Treatment | Maximum new native programs |
|---|---|---|---:|
| Original four development centers | `{water,aquo8}_central_M` | Reuse selected E/analytic forces; retain four full Cartesian checks | 4 |
| Additional four radial development cells | `water_{minus,plus}_M` | Reuse selected E/analytic forces with exact raw-output replay | 0 |
| Consumed hydrated donor transfer | `hydrated_F_S_M` (18) | Analyze plus testgrad each; full Cartesian checks for six central cells, analytic-only for twelve side cells | 36 |
| Unchanged donor controls | `matched_dry_F_central_M` (6) | Prove absence of water-O recipients and reuse exact original E/forces | 0 |

**Maximum: 40 new native program invocations, ten Cartesian checks and 1,380 internal displaced energy evaluations.** The four additional water cells are now development data, so their old radial-transfer label must not be interpreted as transfer or independent validation. The original four-center subset remains separately reported for continuity.

This preserves the previous ten numerical derivative checks rather than adding a new derivative panel. The new four radial development cells have exact analytic-force replay; no new Cartesian qualification is claimed there. All eight water rows must be admitted before the donor stage. Eligible prior complete results with identical physical inputs, engine, settings and raw outputs reduce the maximum; never repeat an identical completed calculation. Reusing a prior full Cartesian result requires the same full numerical evidence, not merely an analytic gradient.

Before donor work, all four original-center Cartesian gates must pass: pinned native `testgrad system.xyz Y Y 0.00001`, total separation 1e-5 Å (individual displacements ±5e-6 Å), every component present, and maximum analytic-minus-numerical gradient component ≤0.001 kcal/mol/Å. Center E/force replay tolerances remain 0.0001 kcal/mol and 0.001 kcal/mol/Å. Apply the same Cartesian gate to the six hydrated centers. Preserve the finite/native/polarization gates, first-failure scheduling stop, all sibling receipts and full32 denominator. No retry or changed derivative step follows failure.

## Minimal adapter changes

Use a new sibling adapter and immutable outputs; do not edit the completed17578 adapter or its results.

1. **Reuse existing data and arithmetic.** Load the pinned prepared32 `DATA.json`; keep original32 baseline and reference rows intact. Reuse the existing native parser, rigid-body helper and SCORE arithmetic. The water-only key renderer must still change only represented metal–water-O VDWPR/POLPAIR records. No donor, water-H, charge, polarizability, bonded or permanent-multipole parameter changes.
2. **Change candidate admission explicitly.** Replace the old joint-four-center result/source assumption with the pinned radial source, exact eight-cell selected model and radial objective. Replay that objective using the radial training result, not the original32 baseline as the radial-start baseline. Check the original four-center diagnostic separately. Keep original family-bound origins separate from the new start when reporting parameter bounds.
3. **Follow actual native authority.** A selected radial row may come from a new complete radial cache cell or exact17578 baseline reuse. The latter is not reuse of the original SI baseline. Verify selected effective XYZ/key/prm hashes, saved source identity, executable/helper pins, energy/forces and raw outputs through the correct authority branch. Do not relabel old cache objects under a new evaluator identity.
4. **Reuse all eight selected water rows.** New non-FD water-side calls are unnecessary. Reuse any prior17578 donor result only when its exact effective inputs and native authority match; otherwise perform the declared donor calls. A previous candidate's merely similar parameters are not cache equality.
5. **Keep scoring, update reporting labels.** Preserve every original score record and gate. Add explicit memberships `original_center_development_4`, `radial_water_development_8`, `hydrated_donor_transfer_18`, and `unchanged_dry_controls_6`. Old `water_radial` score labels remain a historical partition only. Report the original SI baseline, previous17576/17578 candidate, frozen radial candidate and fixed PBE0-D4 reference separately; do not replace baseline data in place.
6. **Narrow finite execution accounting.** Set the fresh run's program ceiling to40 and retain internal FD-energy ceiling1380, at most8 CPU workers and one thread per native process. Fresh ledger, protected parent pins, locked result persistence and aggregate child CPU accounting remain unchanged. Saved/native cache reuse is explicitly recorded; an absent reference or failure is never a zero.

Reuse the prior49 saved parser/arithmetic controls without broadening their scientific scope. The necessary adapter checks are candidate-authority branch admission, exact eight-row reuse and score replay, new40-program inventory, unchanged ten-FD inventory, and distinct development/transfer labels. Only these changed interfaces need additional saved checks. No new scientific tolerances are introduced.

## Frozen comparisons and interpretation

All previously declared development resolutions remain unchanged: absolute energy errors ≤1 kcal/mol, each individual and paired per-fixture rigid RMS ≤2 kcal/mol/Å and maximum body component ≤5 kcal/mol/Å. Keep equal rigid-body categories, water averaging within its category, and geometric-centroid torque divided by1 Å. No averaging away a failed metal, family or radial position.

Retain every previous comparison:

- Individual radial works E(side)−E(center), their Dy−La differences and errors. The original water side path moves the metal only outward by±0.15 Å, with derivative −F_metal·u; hydrated side paths move the entire ligand, with derivative −sum(F_ligand)·u. The0.30 Å secant remains a coarse shape diagnostic.
- Balanced exchange `[E_Dy(q)−E_La(q)]−[E_Dy(aquo8)−E_La(aquo8)]`, using each model's own aquo8 denominator.
- Central matched water addition `[E_Dy(wet)−E_La(wet)]−[E_Dy(dry)−E_La(dry)]`. Free-water and element offsets cancel across the paired metals. This is a metal-dependent many-body water-addition contrast, not a unique pair interaction.
- Individual La, Dy and paired rigid-force errors, with internal atomic forces retained separately as diagnostics. Report unchanged dry failures independently from water-induced changes in hydrated donor behavior.

Use the same uncorrected PBE0-D4 energies and gradients. CP/CC uncertainty remains a separate method diagnostic, never mixed into the objective or substituted for a failed check. Preserve original states: charges+3 for water/aquo8/amide and+2 for acetate, physical La1/Dy6 labels and saved ECP effective states.

| Outcome | Interpretation |
|---|---|
| Any native/Cartesian gate fails | Numerical or admission failure; no transfer qualification |
| Eight water cells pass; hydrated donor cells fail | Water development succeeds; coupled donor/water mechanics remain inadequate |
| Hydrated errors partly improve but some frozen checks fail | Partial consumed-data transfer with remaining failures; no promotion |
| Unchanged dry controls fail | Existing donor-model inadequacy; not deterioration caused by water parameters |
| Every numerical and consumed-data check passes | Useful development transfer on these configurations only; no independent validation, bulk-water, protein, periodic or affinity claim |

The eight water cells have all influenced parameter selection. They cannot be used to claim improved generalization after this radial refit. The donor configurations are consumed development-transfer data; the reserved24 outcomes remain unopened.

## Existing pins

Paths are relative to the project root; full paths and hashes are also in `PINS.json`.

- Previous transfer specification/source/data/score/tests: `diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_transfer_v1/`.
- Prepared32 inventory: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_transfer_v1/prepared_v1/DATA.json`, SHA `1f8e85d3349c447b61ec8f40b9ec8be56c9ecab02f0926272180d16237d4af10`.
- Previous complete transfer: `.../paired_water_transfer_v1/run_v1/execution/RESULT.json`, SHA `702d97ababc232336e3b5405a1326cec4cafeccc9570ce9f514076cf7908e986`.
- Frozen radial SPEC: SHA `0e937e4ef9ddf3d768ff93509c5e8e705a5384956895fa3e46efc89e8bda7f41`; source REFINE `e808ae45f65380255cb02dbace2865932972873989ab3ee8782ccb5f050b22e0`; objective `5b83cf5f3e77aeb96ec0f4d0a0c3f765883f406e8f4110648972e5bac0996f6b`; evaluator `1cf6025e27291523ac7c7242a07245da22c8927d67f4323c87392974133eda56`.
- Radial training: `.../paired_water_radial_repair_v1/training_v1/RESULT.json`, SHA `3aefeb74f40aaf51e327fdc5c0ffd85f7f00e6e769b07c080f3ece16c345192d`.

Candidate result, bound inputs, runnable config and root release: **unassigned**. This plan launches nothing.
