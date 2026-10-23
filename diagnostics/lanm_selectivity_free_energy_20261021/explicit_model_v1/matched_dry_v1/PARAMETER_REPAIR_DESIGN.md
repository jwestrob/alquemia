# Finite paired-ion repair proposal

Proposal only: no fitting, new reference calculation, geometry scoring or dynamics is released here. First finish the six matched-water-removal cells and their saved comparison. Existing water, dry-dimer, hydrated radial and matched-dry panels are **consumed development evidence**, including failed or unfavorable observations.

## What to repair and what to keep separate

The individual force error is largely shared ligand internal distortion; fitting ion parameters to that total RMSE would target the wrong observable. Fit balanced exchange energies and metal/rigid-molecule forces and torques. Retain full atomic forces as diagnostics, including the paired internal response, without pretending the reference uniquely partitions bonded and nonbonded physics.

The current Dy ion was transferred from a different water parameterization, while La has published donor-specific corrections. A donor fit must not silently compensate for an incompatible ion–water model. Keep charges, intrinsic ion polarizabilities, ligand/protein permanent multipoles, bonded terms and water geometry model fixed. No fit to LBT preference, experiment, optimized protein structures or the inspected desired sign.

## One small staged parameter family

Use the native `VDWPR` and `POLPAIR` mechanisms, with exact atom-type/class mappings and analytic forces. La's published model remains a fixed control; its absolute disagreements remain visible. Failure of that control cannot be cured by adjusting Dy to match a difference.

1. **Dy–water stage: three parameters.** A water-O/Dy pair radius, pair well depth and pair Thole damping coefficient. These pair-specific overrides leave the intrinsic ion parameters and every donor pair unchanged. Train only against admitted existing water-dimer and aquo8 paired data. Include the balanced dimer-versus-aquo8 exchange and rigid-water/metal force response, not the shared 104.52° water bending residual.
2. **Dy–donor stage: four additional parameters.** One pair radius and one pair Thole coefficient for carboxylate oxygen, and the same two for amide oxygen. Freeze the water-stage result and all pair well depths at their inherited values. The same carboxylate pair must serve mono- and bidentate acetate; no geometry-specific corrections. Train against all admitted central dry/wet donor pairs, balanced aquo-referenced contrasts, water-addition double differences and rigid ligand forces/torques. Radial development data with no qualified reference remain diagnostics rather than invented fitting targets.

This is a deliberately restricted seven-parameter maximum. The donor well-depth restriction may prove inadequate; reject that family rather than adding parameters after viewing validation. Before a fitter runs, bind numerical bounds, energy/force weights, optimizer budget and reference-accuracy requirements in one executable specification. A reasonable starting trust region is ±15% of the inherited pair radius, water well depth 0.25–4 times its inherited combined value, and damping 0.05–0.50. These are proposed constraints, not established physical intervals. Repeated bound hits, weak parameter identification or loss of the water control are stop conditions, not grounds for automatic expansion.

Only reference cells with the declared electronic-state and derivative credibility may become fitting targets. Successful native termination alone is insufficient. The current frozen-f PBE0 comparisons still carry that qualification boundary.

## Minimal new validation matrix to freeze before any fit

Prepare and hash these coordinates before optimizing parameters; keep their reference/model outcomes unopened until candidate selection:

| Validation family | New geometries | Both-metal cells |
| --- | ---: | ---: |
| Aquo8: two rigid-water orientational perturbations, oxygen positions fixed | 2 | 4 |
| Hydrated mono-acetate, bi-acetate and acetamide: two rigid-ligand orientations per family, waters fixed | 6 | 12 |
| Existing prepared ACE–Asp3/Glu11/Asn5/Trp9 protein donors | 4 | 8 |
| **Total** | **12** | **24** |

For the first eight new configurations, use opposite ±15° rotations about deterministic geometry-defined axes, without changing covalent geometry. Freeze the exact axes/centers and contact gates in preparation; a failed contact gate remains an unavailable validation slot, not a selectively replaced geometry. For bidentate acetate use the coordinating-oxygen midpoint as center; use the coordinating oxygen for mono-acetate/amide. For aquo8 rotate each water about its oxygen using a defined local tangent axis. Both metals share each coordinate set. No new water or ligand counts, optimizations or radial search are required.

Report every validation cell, absolute and paired errors, and the reserved orientation response; do not fit again after opening them. An improvement confined to development energies, a force-sign regression, or compensation between individually poor La/Dy responses does not qualify the model. Concrete numerical release thresholds must be declared with reference uncertainty before fitting; the existing 1 kcal/mol reporting resolution is not a reference error bar.

## Necessary protein transfer test

The eight capped protein geometries are already prepared and unscored, with exact native protein donor multipoles/polarizabilities. They are a transfer test, not another calibration set. Make the complete mapping explicit: SI carboxylate O type/class403 maps to protein Asp143/Glu161 (class31); SI amide O type105/class55 maps to protein carbonyl class5, including Asn153, backbone5/11 and ACE cap227. SI amide carbon type104/class53 maps to protein carbonyl class3, including backbone9, Asn152 and ACE226. The published La carbon `VDWPR` must be carried through explicitly alongside its oxygen correction; omitting it would change the baseline.

Use classes for `VDWPR` and explicit chemical-family type lists for `POLPAIR`, including all represented cap carbonyls rather than selecting only a favorite coordinating atom. Freeze that complete map before scoring. Protein donor multipoles/polarizabilities remain unchanged and differ materially from SI values; semantic mapping alone does not qualify the transferred correction.

Apply a selected donor correction to these protein types only as a declared hypothesis. Score all eight with the paired reference convention and report metal forces, ligand rigid response and full atomic forces. If transfer fails, retain it as a failed validation and require a new independently validated protein-donor parameterization. Do not fit the protein fragments retroactively or proceed to production sampling. The existing whole amidated LBT3 box remains a topology/chemistry preparation, not a validated Hamiltonian.
