# Archived donor-response probe, frozen before feature extraction

Authorized parallel structural-response experiment, resumed 2026-09-19 after
Jacob's goal discussion. Parent approved these exact arms and archive-only scope.
Protocol: `mace_polar_archived_donor_response_probe_v1`. Research only; no default
or historical preparation/result changes, no new molecular evaluations.

## Question and population

Do Ca-minus-La donor gradient patterns provide class information alone or add to
the existing native DFT contrast? Reuse 56 native MACE-POLAR-medium core endpoint
forces: all 25 canonical PQQ controls and the three consumed crystal transfers
1H4I, 4MAE, 1KB0. Separately reuse eight native MACE and eight native DFT analytic
gradients for GGR extended/connected (both 1GLG) and alpha-lactalbumin 1F6S/6IP9.
The direct population has two biological groups, not four independent proteins.
All cases are development/retrospective, never prospective blind validation.

Pinned archive roots: `workspaces/mace_canonical_20260916/mace_v2/medium`,
`workspaces/mace_canonical_20260916/audit_v2/inventory.json`,
`workspaces/mace_metal_response_20260918/prepared_v2/preparation.json`, and
`workspaces/site_classifier_20260918/{features_v2,evaluation_v1}`.

## Fixed physical observables

Delta g = g_Ca - g_La, with native MACE g = -F times 23.06054783061903
kcal/mol/eV. Native DFT g uses the released Hartree conversion and bohr length.
Each pair has identical nuclear coordinates, atom order, microstate and water
inventory; charge differs by one. The primary MACE Hamiltonian is the identical
archived float64 native MACE-POLAR-medium vacuum Hamiltonian across cases. It is
not the masked OMOL model, CPCM, or the MACE+GB response model. DFT is the separate
native r2SCAN-3c/CPCM Hamiltonian: comparisons do not establish force equivalence.

Apply existing chemical donor typing to retained physical source atoms, then a
common 3.2 Angstrom metal distance cutoff. Use exact schema-specific PQQ allowed
donors, backbone/sidechain donors and existing explicit water oxygens. Canonical
archived donor lists are truncated at 3.1 Angstrom and cannot implement this
selection directly. No new atoms, donors or waters are added to the core.

For e_i = (x_i-x_metal)/distance, the two primary features are:

1. Mean radial differential gradient: mean_i(Delta g_i dot e_i), kcal/mol/A.
2. Tangential fraction: sum_i |Delta g_i - e_i(Delta g_i dot e_i)|^2 divided by
   sum_i |Delta g_i|^2. Undefined zero-load cases remain unavailable.

Donor coordinates are physical Cartesian displacements with other physical
atoms fixed. Synthetic link hydrogens are eliminated through their exact cap
chain rule onto retained AND omitted source atoms. Real physical hydrogens stay
explicit; they are not synthetic mechanical modes. Donor-load coherence
|sum Delta g_i|/sum|Delta g_i| and complete-PQQ net differential force and torque
about its geometric centroid are descriptive only, never extra fitted features.
These do not estimate equilibrium structure, scaffold stiffness or uncertainty.

## Fixed evaluation

Use the existing four canonical homology folds (18,5,1,1), class/group-balanced
weights, training-only standardization, logistic ridge lambda=1 and zero-logit
class rule. Five fixed arms, all reported: DFT_R only; radial only; DFT_R+radial;
tangential only; DFT_R+tangential. At most two features; no tuning or feature/arm
selection. R differs from published S by one common aquo offset in this stratum,
so standardized DFT-only predictions must reproduce the earlier S-only model.
Compare to archived DFT+structure (25/25) and DFT+structure+masked-MACE (24/25).
No added-accuracy claim if the existing DFT-only arm already calls all 25.
Fit all 25 only for the three consumed crystal checks; report known homology
overlap and leave 1KB0's uncomputed validation group explicitly unknown.

For direct cases report features/orderings, donor-vector cosine agreement with
DFT, and GGR partition differences. No affinity class from arbitrary feature
signs, no PQQ threshold transfer, no classifier fit on the two direct groups.

## Essential checks and execution

Verify input hashes, exact paired XYZ/charge/order, native model/software match,
force shapes, analytic DFT input/output/energy/coordinate identity, source maps,
and complete PQQ inventory. Check observable rigid-rotation/translation algebra
(1e-9 feature tolerance), cap-chain-rule force/torque closure (1e-4 kcal/mol/A or
kcal/mol, accommodating six-decimal cap coordinates), and source-heavy coordinate
agreement (1e-6 A). These are archive/parser/algebra checks, not new scientific
integration evaluations. Failures remain explicit with their original denominators.
Use one CPU process, OPENBLAS_NUM_THREADS=1, existing lanm_qmmm Python. Record
wall/CPU time and peak RSS. No new GPU, DFT, optimization or molecular-model calls.

