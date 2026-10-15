# Asn290 amide-H provenance and already completed checks

## Exact current source

The reported1.183–1.188Å bonds are **Asn290 ND2–HD21/HD22**, not backbone N–H. These are modeled physical protein hydrogens added during preparation, not raw AF3 hydrogens, PQQ-carver hydrogens or synthetic boundary caps. Ca4's raw CIF and both Ca4/La4 normalized PDBs contain only the eight Asn heavy atoms; HD21/HD22 first appear in the corresponding protonated PDBs.

The September20 source preparation is recorded in:

- [Ca4 protonation manifest](../../../workspaces/accommodation_goal_20260920/folds_v1/source_preparation/cases/a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4/protonated_protonation_manifest.json), SHA256 `f1215683d6fd66ff835ce9a53f9b24df86125745bb90038dd9604dc2666cfec0`.
- [La4 protonation manifest](../../../workspaces/accommodation_goal_20260920/folds_v1/source_preparation/cases/a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-4/protonated_protonation_manifest.json), SHA256 `0a75b50b4c6a778593ae46cff1d107638a939d00d6bf7556ed38d2e1d4017fd6`.

Both use `pdbfixer_standard_only_rng20260914_openmm_cpu_threads1_v1`: PDBFixer1.12.0/OpenMM8.5.1, pH7, RNG20260914, CPU1, `forcefield=None`; no missing heavy atoms or residues were rebuilt. The pinned wrapper `protonate_standard_only.py` SHA5c0393c1… and canonical `scripts/protonate_cif.py` SHA9a1adf14… still match. Current204-atom maps inherit those source coordinates. The amide flip and subsequent torsional searches preserve these bond lengths.

## What “minimized” means here

PDBFixer calls `Modeller.addHydrogens(forcefield=None)`. Matching-version source builds H-parent harmonic restraints plus generic repulsion, fixes existing atoms, and calls `LocalEnergyMinimizer.minimize(context,1.0,50)`. Thus these H coordinates underwent the builder's generic minimization; no exact-source iteration/force-convergence receipt was captured. They were not minimized with Amber, MACE or DFT during this preparation.

This mechanism was **already measured**, not merely inferred from today's lengths. The [September20 diagnostic](../../accommodation_controls_20260920/HYDROGEN_DIAGNOSTIC_REPORT.md) captured Systems with zero nonbonded exclusions, including bonded pairs. Its recorded isolated H-parent potential minimum is1.17705Å. Successful MMOL1770 La3 and repaired La4 had median N–H1.185/1.186Å. Continuation1203809 reached its declared force criterion while retaining1.170993–1.207169Å parent-H lengths. That repaired collapse/convergence under the **same** potential; it did not establish physical bond lengths. Exact A0AC H-length correction was not found in the reviewed archives.

## Existing alternatives: what was actually validated

| Archived route | Actual result and relevant limit |
|---|---|
|[Protein-H radial projection](../../mace_hydrogen_20260916/REPORT.md), `scripts/mace_hydrogen.py::repair/validate`|Uses actual ff19SB topology equilibrium lengths for every protein C–H/N–H/O–H, preserves directions/heavy/cofactor coordinates, and verifies prescribed lengths and identity on real1H4I fixtures. Four full-protein MACE calls did not resolve model disagreement. This bypasses generic placement; it does not modify its bonded repulsion or validate an A0AC correction.|
|[PQQ contextual physical-H optimizer](../../second_shell_20260919/PQQ_H_REPORT.md), archived `pqq_h_v1/implementation/second_shell_hydrogen.py`|Native OMOL optimizes physical H with fixed heavy atoms/caps and parent/displacement gates. Four1H4I/4MAE proposals and four DFT checks lowered endpoint energies; DFT separation30.571889→30.071899kcal/mol gave no accuracy gain. These are different sources and constrained preparations, not fixed bond-length targets.|
|[Contextual water](../../contextual_water_20260919/REPORT.md) / [coordination](../../coordination_preparation_20260919/REPORT.md)|Water-only preparation leaves all25 dry canonical PQQ inputs unchanged. Coordination starts from the crystal H-prepared geometries and tests heavy modes; neither tests A0AC amide-H lengths.|
|[September28 LanM repair](../../metal_environment_response_20260926/scaffold_H_repair_20260928/REPORT.md)|Native bonded+steric optimization improved C-bound H geometry on Hans/Mex; exchangeable N/O-bound H stayed fixed. It supplies no A0AC ND2-H validation.|

No archived variant removing bonded-pair repulsion from the generic OpenMM placement System was identified. Tested alternative coordinate methods exist, but neither has been validated on these exact A0AC states. The concrete provenance warrants retaining this representation limitation; it does not establish that it causes the rotamer response. Job4971 remains the matched inherited-H experiment. No correction, calculation or submission accompanies this note.
