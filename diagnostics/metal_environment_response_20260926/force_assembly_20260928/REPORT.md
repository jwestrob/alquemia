# Finite physical gradient assembly: real archived mapping works

**Archived La A/B electronic gradients now map onto all1,891physical source atoms.** Both cap anchors and all1,696external point-charge gradients are included. Five real-artifact mapping/parser tests pass. No electronic or classical molecular evaluation was performed.

This tests an assembly operation using the **old native r2SCAN-3c, original-H Hans gradients**. It does not substitute their energy/gradient into the repaired frozen-4f experiment. The code rejects mismatched core/field source hashes, atom order, coordinates, energy or gradients. Missing endpoints remain unavailable.

## Expression and units

For cap position c=r+L(o−r)/|o−r|, the physical anchor contributions are J_r^T g_c and J_o^T g_c, with J_o=L/|o−r| (I−uu^T), J_r=I−J_o. The code recomputes and checks both Jacobians and the cap coordinate from actual source anchors before use.

Real QM gradients scatter to their source IDs; the exchanged target uses its real source-site ID (`A/203//LA` for this source even when its endpoint element is Dy). Point-charge gradients scatter in verified serialized field order, including spectators. Synthetic caps never become physical atoms. Boundary omitted CA atoms remain present and receive the omitted-anchor contribution even when their point charge is removed.

Native engrad and pcgrad are dE/dR in Hartree/bohr. Convert exactly once to kcal/mol/angstrom. Forces are the negatives. With a compatible classical gradient, `g_total = g_electronic_mapped + g_classical`. This module computes no Coulomb energy/force, so it adds no second QM–MM electrostatic term. Fixed redistributed charges are treated as constant parameters; no coordinate-responsive charge model is asserted.

The classical adapter requires explicit source IDs, coordinates, gradient units/sign and exclusions of QM–MM Coulomb/C4 induction. It reorders by ID and rejects incomplete inventory or mismatched coordinates. Without that actual matching classical artifact, the combined gradient/force remains **null**, not the electronic component relabeled as a complete hybrid result.

## Executed evidence

Both A/B archived results map195QM centers (including4caps) and1,696point charges to1,891physical atoms. The change from the original sum of native QM+point-charge gradients to the mapped physical sum is at most1.17e−13kcal/mol/angstrom per Cartesian component. Actual translation residuals remain approximately5.5e−6kcal/mol/angstrom in norm and are reported without subtraction. This algebraic conservation does not establish numerical exactness of the original Hamiltonian.

Tests check both real configurations, both-anchor chain rules using actual saved cap gradients contracted with perturbed source mappings, rotation covariance of the mapping algebra, wrong-source/missing-endpoint rejection, and forbidden classical Coulomb flags. The linear contraction used for chain-rule finite differences is explicitly an algebraic functional, not a new molecular energy or synthetic scientific validation dataset.

Compact result inventory: `ARCHIVED_MAPPING_RESULT.json`. Full mapped gradient arrays stay under `workspaces/metal_environment_response_20260926/force_assembly_v1/`.

## Runnable interface

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/metal_environment_force_assembly.py --inputs workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8DQ2/INPUTS.json --configuration A --metal La --collection workspaces/metal_environment_response_20260926/lanm_ef3_hans_scout_v1/FINAL_COLLECTION.json --task-id La_A --output workspaces/metal_environment_response_20260926/force_assembly_v1/ARCHIVED_LA_A_reproduction.json
```

Optional `--classical` accepts the explicitly declared bridge artifact; it does not run mechanics. All output writes are exclusive-create. Real repaired gradients must be supplied with their matching preparation and collection. No affinity, whole-protein relaxation, solvent consistency or validated full hybrid Hamiltonian is claimed. The final ledger's boundary convention and model assumptions remain separate scientific checks.
