# Paired water response: four-parameter representability pilot

## Question and scope

Can separate La and Dy water-contact radii and polarization damping reproduce the consumed uncorrected PBE0 energy/force response without hiding either metal's error? This is a new, bounded **representability experiment**, not a physical calibration or sampling release. The held Dy-only seven-parameter family stays unchanged and unsubmitted. Root selected four parameters before fitting; a six-parameter family with adjustable epsilons is an unrun future hypothesis, not an automatic fallback.

Use exactly `water_central_La`, `water_central_Dy`, `aquo8_central_La`, and `aquo8_central_Dy` from the existing training result. Keep coordinates, bonds, water geometry, multipoles, charges, ion/water polarizabilities, all covalent terms, all non-water pairs and mutual-polarization settings fixed. There is no implicit solvent, geometry search, donor fit, protein score or new reference call. This uses the small-molecule SI water model, including its 1.5337 Å Urey–Bradley reference distance; the protein model's 1.5537 Å value is a separate transfer issue.

The target is the existing **uncorrected PBE0-D4/def2-TZVP frozen-f ECP energy and gradient surface**, La core46/Dy core55, effective singlets and physical labels La1/Dy6. Re-admit exact native/state/geometry/source evidence and complete finite arrays for all four cells. Water48 supports the sampled dimer directional derivative and shows the La mismatch survives the tested method/CP choices; it does not qualify every aquo8 force, explicit-f physics, bulk hydration or the reference as truth. Published La calibration used CP-DLPNO-CCSD(T)/ZORA curves and hydration free energy; this pilot can move away from that calibration. No CP/CC energy is mixed with uncorrected PBE0 forces in the objective.

## Four parameters, unchanged context

Parameter order is `(La_R, La_damping, Dy_R, Dy_damping)`. For each metal, change only water-O `VDWPR` radius (class34 paired with class500 or501) and water-O `POLPAIR` damping (type36 paired with type500 or501). Retain each pair's resolved initial epsilon exactly. Water-H pairs and the atom-level ion parameters remain unchanged; there are no new polarization-group bonds.

Resolve initial R, epsilon and damping using the pinned `alchemical_path_v1/PATH.py::parse_parameters/resolve_pair` against the exact original parameter/key text. This accounts for explicit overrides and native CUBIC-MEAN/DIAMETER and HHG conventions; atomic diameters are not pair radii. La's initial resolved damping is **0.25**, Dy's **0.139**. Radius bounds are `[0.85 R0, 1.15 R0]` independently for each metal; both damping bounds are `[0.05, 0.50]`. Bounds are exploration limits, not uncertainty estimates, and cannot expand during this run. Preparation must publish the resolved numerical family and prove the initial vector resolves to the original Hamiltonian before baseline reuse.

## Symmetric objective and separate controls

Let `D(q)=E_Dy(q)−E_La(q)` and `C=D(water_central)−D(aquo8_central)`. C is the balanced reaction energy for `DyW8 + LaW1 → DyW1 + LaW8`; both sides contain one ion of each element and nine waters. Use only its model-minus-reference error, scaled by 1 kcal/mol. Do not fit raw unlike-metal total-energy offsets or unbalanced dimer-to-aquo8 energies.

For each covalent body (ion or complete water), use its net force and, for water, torque about its geometric centroid divided by 1 Å. Average squared Cartesian residuals, weight translation/rotation equally within each water, average waters within their category, and weight the ion and water categories equally. Average the two geometries equally. Repeat for La, Dy, and matched Dy-minus-La force differences. Forces and torque/1 Å use the same 2 kcal/mol/Å scale. Keep all atomic/internal-force residuals as diagnostics; do not fit stiff internal modes through ion parameters.

The scalar loss is **0.25 C-error² + 0.25 La rigid-response loss + 0.25 Dy rigid-response loss + 0.25 paired rigid-response loss**, using those scales and normalized block means. The individual and paired blocks share information; these are explicit development weights, not independent likelihood observations or statistical confidence.

Report separate development checks: `|C error| ≤1 kcal/mol`; each of La, Dy and paired category-balanced rigid RMS `≤2 kcal/mol/Å` and maximum body Cartesian-component error `≤5 kcal/mol/Å`. Require all checks to label a candidate development-consistent. A smaller scalar loss alone does not pass them. Both metals' results remain visible; neither can be excused by cancellation.

For baseline and selected candidate, project the existing dimer forces along the exact Water48 outward-water direction and report individual/paired residuals against all four saved QM derivatives (uncorrected/CP PBE0-D4 and uncorrected/CP CC). Flag incorrect derivative signs and any absolute CC residual that worsens relative to baseline, separately for each metal and the pair. These sensitivity flags do not alter the search or substitute CC targets; retain the full magnitudes. They preclude claiming the candidate repairs reference-method uncertainty. No candidate receives physical qualification from this pilot, even if all development checks pass.

## Finite search and identifiability

Use deterministic coordinate pattern search in normalized intervals: initial step0.10; poll minus then plus for each of the four parameters; select the lowest objective only after the complete poll, with fixed-order ties retaining the earlier point; halve after a complete poll with no improvement; local stop when the next step is below0.00625. Clip bounds and deduplicate exact repeated vectors. Parallel completion order cannot affect selection.

Allow **241 distinct vectors total**, including the reused baseline and **eight reserved final sensitivity vectors** (±0.02 normalized interval units per parameter, clipped at bounds). Thus the search can use233 vectors including its baseline, enough for29 full eight-point polls before clipping/reuse. Do not start a poll that exceeds the remaining search allowance; return `budget_limited` with the best completed point. Budget exhaustion does not establish failure of the parameter family or optimizer convergence. Preserve the reusable exact cache.

Sensitivity vectors cannot replace the selected candidate. Report weighted-residual Jacobian column norms, relative1e−6 rank, condition estimate and flag condition>1e6 or rank<4; report proximity within1% of a parameter bound. These are local identification diagnostics, not parameter uncertainty or a global optimum. A flagged candidate cannot be called uniquely identified. No automatic family/bound expansion follows.

## Execution and release boundary

Implementation may reuse the pinned parser, body-response helpers and native `analyze E`/analytic `testgrad Y N` admission. It must independently compare the original four baseline records and publish exact geometry, key, parameter, executable and source hashes. Full-vector request identity includes all four parameter values. Native-cell cache identity instead includes the cell, its rendered physical input and complete model/source lineage: a La-only parameter change reuses unchanged Dy cells, and vice versa. Preserve each request's association with its native-cell authority, including proven baseline equivalence. No completed identical native input is repeated. Conservative ceiling: **1,920 new native program invocations** (240 nonbaseline vectors ×4 cells ×2 programs); cache reuse can reduce this. Preserve four-cell availability and all attempted receipts on failure.

Root owns one8-CPU `standard_expand` allocation. Run at most8 independent complete-poll vector workers, each with one native process at a time and all numerical thread controls1; no login fitting or GPU is needed. First native failure/nonfinite/polarization admission failure stops new scheduling; retain already-running siblings and their results, with no solver fallback. Record actual program, CPU and wall costs. Candidate force finite differences, hydrated development transfer, reserved24 validation and periodic transfer are later separately declared gates, not hidden calls in this budget.

This file freezes the scientific design for source implementation/review only. Root binds an actual execution release after source and saved-data admission review. No fitting, reference calls, validation access or dynamics has occurred in preparing this specification.

## Frozen evidence

- Four rows within training result: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/matched_dry_v1/training_v1/RESULT.json`, SHA256 `99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb`.
- Water48 collection: `workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/water_force_qualification_v1/collection_v1/RESULT.json`, SHA256 `7be6763c6eddaa1501c0e895851433ab7117b6e039472df5b0bf50561a90e818`.
- Pair resolver: `diagnostics/lanm_selectivity_free_energy_20261021/alchemical_path_v1/PATH.py`, SHA256 `69e1e4dc46cce86f145596180a4507d1e09b3e6473a83210efd99d19e16abba3`.
- Original key SHA256 `6cdcc95d348953b8aafdb11b89cc03946eafee20eaa9a301622f8a8c3715cd89`; original small-molecule parameter file SHA256 `e1cc2bcc6d7fb319aba7abaaf6b4483f2ef5ba730727212b5c157a18a290f38b` (paths in the pinned training lineage and prior family).
- Periodic mapping and water-model distinction: `diagnostics/lanm_selectivity_free_energy_20261021/periodic_mapping_v1/REPORT.md`, SHA256 `c807826445664953e0337c2bdeab115b644ccfe14b2906202b59b67932a1d1c5`.
