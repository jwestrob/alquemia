# Hans 8DQ2 EF3 finite-field preparation

Ready for the separately owned four-endpoint La/Dy electronic-response scout. No molecular energy, force, optimization, scheduler submission, or library outcome access occurred here. The exact real preparation passes seven mapping/state/geometry checks. This does not qualify the electronic surface or a complete QM/MM model.

## Frozen source and question

Use the normalized whole-chain source already audited in `../lanm_preparation_feasibility/RESULT_v2.json`, not either failed relaxed candidate. Chain A is a conditional monomer, not an assembly or folding equilibrium. The region is the preselected complete EF3 loop 83–94, the adjoining C/O82 and N/H95 peptide units, supporting C/O65–N/H66 peptide unit, and full waters337/350. There are190 real nonmetal atoms, one exchanged metal, and four caps:195 QM atoms. The physical source retains all1,887 protein/water atoms and all four observed ion sites. Seventy-four crystal waters remain; two are QM.

A is the normalized source. B applies the already declared right-handed +2° rotation about Asp85 CA→CB to CG,OD1,OD2,HB2,HB3, preserving CB and all other atoms. Both targets receive both exact geometries, charge−1, and the identical finite external field. LaIII uses multiplicity1; DyIII uses the physical high-spin sextet hypothesis, multiplicity6. All-electron counts are852/861. Backend ECP/core-electron and open-shell qualification belongs to the reference runner; no effective-singlet substitution occurs here.

EF1/EF2 remain LaIII monopoles+3 each; EF4 remains the observed Na+1, not a fourth lanthanide. These frozen spectator charges are an explicit approximation, with no ion force-field or polarizability claim.

## Native field and local boundary model

The exact1,887 source IDs map to the original protonated topology and native OpenMM8.5.1 ff19SB/TIP3P parameters. No hydrogen addition, atom repair, ion typing, Context, or energy evaluation is used. Charges attach to the archived normalized coordinates. Parent protein/water charge is−4.

All selected QM charges and each cut's MM1 CA charge are removed from the field. Those four CA nuclei stay present in the physical mapping with zero field charge; they anchor the cap mapping. For each partial residue, its remaining MM atoms retain the original full-residue force-field charge because the selected capped peptide fragment is formally neutral. The necessary difference is split equally over the two actual bonded exterior heavy neighbors of the omitted CA:

| Partial residue | Removed atoms | Recipients | Increment per recipient/e |
|---|---|---|---:|
| Leu65 | C,O,CA | N,CB |−0.01120|
| Glu66 | N,H,CA | C,CB |−0.09150|
| Ala82 | C,O,CA | N,CB |+0.03155|
| Trp95 | N,H,CA | C,CB |−0.08565|

This **new explicit residue-local heavy-neighbor rule** is graph checked. It conserves local charge, not the molecular dipole or original electrostatic potential. Before/after MM dipoles are recorded as representation diagnostics. It is not asserted equivalent to native ORCA charge shifting or the previous PQQ N/C redistribution. The full loop contributes formal−4, target+3, giving QM−1; remaining protein/waters sum0 and spectators+7 give field+7. Total physical charge+6 is preserved without global neutralization. The field has1,696 rows; four additional physical MM boundary atoms carry zero field charge.

Caps follow the source cut vectors: C–H1.09Å and N–H1.01Å, with analytic derivatives to both actual retained/omitted atoms. This fixed-length convention matches the project's peptide-amide preparation, not native ORCA's default ratio-based placement. Caps are not independent physical atoms or mechanical coordinates.

## Actual checks and limitations

`CHECK_RESULT.json` records seven passed real-input checks. Maximum covalent bond-length change is7.8e−15Å; maximum moved-atom displacement.0829489Å. The closest moved nonbonded contact excluding1–2/1–3 neighbors increases2.06896→2.12350Å. All four cap Jacobians agree with coordinate finite differences within6.6e−10. La/Dy coordinates, ordered source maps and state parity agree exactly; A/B field bytes are identical.

Only embedded electronic energy is prepared. Classical internal/cross-region bonded, repulsion/dispersion and solvent total-energy terms are not provided; their absence is not a zero correction. No absolute affinity, source preference, biological classification, force qualification or full-protein relaxation is established. Dy spin-orbit/localization and reference-method sensitivity remain unresolved. Other Hans/Mex sources are not prepared in this delivery.

An initial export-only attempt stopped at an assertion because the source water-kind label is `retained_crystal_water`; its files are preserved at `Hans8DQ2_export_only_attempt1`. The corrected final build made no scientific change to source coordinates or parameters. Final inputs are immutable.

## Files and reproducible preparation

Final contract: `workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8DQ2/INPUTS.json`, SHA256 `c5886704ee541752ec2adfe7ada0b5f8b2c0df6014d5ed8676b74447d1141c3f`.

It pins source artifacts, exact force-field XML and OpenMM export, all four XYZs, A/B core maps, cap Jacobians, local charge ledger, field and environment identities. `source_state` points to the actual pre-scoring feasibility inventory.

From repository root, an unused output directory is required:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/lanm_ef3_preparation/PREPARE.py --repository "$PWD" --output workspaces/metal_environment_response_20260926/lanm_ef3_preparation_reproduction/Hans8DQ2
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/lanm_ef3_preparation/CHECK.py --repository "$PWD" --inputs workspaces/metal_environment_response_20260926/lanm_ef3_preparation_reproduction/Hans8DQ2/INPUTS.json --output workspaces/metal_environment_response_20260926/lanm_ef3_preparation_reproduction/CHECK.json
```

Root owns the subsequent `metal_environment_lanm_reference.py` manifest, state support, scientific calls and collection. Preparation used no GPU or scheduled compute; no unmeasured CPU utilization or speedup is claimed.
