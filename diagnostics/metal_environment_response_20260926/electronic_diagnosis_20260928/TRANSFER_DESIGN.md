# La/Dy EF3 transfer design — 28 September 2026

## Decision and observable

After the current frozen-f engine/derivative gate, a useful finite comparison is **both consumed Hans source states against the same consumed Mex source**, with no source selection or averaged-away reversal. There is one protein-level qualitative comparison: Hans has stronger relative La-over-Dy apparent affinity than Mex. There are no independently measured EF3 affinity labels.

For source/state `H` and `M`, define

`D(H,M) = [E_H(Dy) − E_H(La)] − [E_M(Dy) − E_M(La)]`.

This is the energy of `H–La + M–Dy → H–Dy + M–La`. Positive D supports stronger *relative* La preference in Hans on this conditional electronic model. Every protein, water, proton, spectator and cap occurs once on both sides: the two proteins need not have identical formulas for this exchange to balance. Absolute totals of unlike regions are nevertheless never alternative geometries of one system. No raw cross-source minimum is allowed.

Use one explicitly declared Hamiltonian, matching per-element basis/ECP assets, solvent convention and energy expression for every cell. Frozen-f Dy ECP55 is a different model from native explicit-f r2SCAN-3c: no old Dy endpoint, native 3c offset or MACE/GFN2 total can fill a new cell. La needs a deliberately specified compatible treatment before this becomes runnable discrimination. The atomic element/core offsets cancel within the balanced cycle only when treatments are fixed across both proteins. Physical DyIII remains4f9/S=5/2; effective restricted valence multiplicity1 must stay separate from that physical state.

## Actual source compatibility

The following counts come directly from existing endpoint XYZ and pinned `INPUTS.json`, not from reconstructed prose. Formulas include synthetic caps and the target Ln.

| Representation/source | Formula | Atoms | QM charge | QM waters | Spectators represented |
|---|---|---:|---:|---:|---|
| Older compact Hans8DQ2 EF3 |C13H24N2O10Ln|50|−1|0|none|
| Older compact Hans8FNR EF3 |C13H24N2O10Ln|50|−1|0|none|
| Older compact Mex8FNS EF3 |C10H21NO11Ln|44|−1|2|none|
| New embedded Hans8DQ2 EF3 |C56H91N17O29SLn|195|−1|2|EF1LaIII/EF2LaIII/EF4NaI|
| New embedded Hans8FNR EF3 |C56H93N17O30SLn|198|−1|3|EF1DyIII/EF2DyIII/EF4DyIII|
| New embedded Mex8FNS EF3 |C56H102N17O33Ln|209|0|7|EF1NdIII/EF2NdIII/EF4NdIII|

The two older Hans cores share composition, donor/cap identities and no explicit waters; their mapped geometries differ. A within-metal source energy difference is composition-balanced there. The compact Mex core is not composition-compatible for pooling with either Hans core, but the four-cell protein–protein exchange above remains balanced. These compact cores exclude the exterior; they test local chemistry, not spectator effects or scaffold relaxation.

The new Hans pair differs by one QM water, total crystal water inventory74 versus179, spectator identity/charge, and finite field. Their apparent geometric transfer also changes chemical context: **do not pool195/198 atom endpoints or call their raw energy difference deformation work**. Mex has167 source waters, a distinct sequence/core formula and charge. All three embedded pairs are individually fixed-composition La/Dy substitutions, so separate cross-protein exchanges remain defined but conditional on their distinct source environments. No water exchange occurs inside that exchange cycle; comparing alternative occupancies as states of one protein would require reservoir accounting separately.

New source maps retain complete loops Hans83–94/Mex84–95, complete adjoining peptide units and nonlocal amide support Hans65–66/Mex66–67. Four source-connected caps are retained. A/B perturbations are the predeclared +2° whole Asp85/Hans or Asp86/Mex sidechain rotation, including HB2/HB3. A/B have identical fields within each source. Source-dependent water and spectator differences are not corrected by deleting atoms or relabelling EF4.

All structures use actual author chain A, conditional monomers. They are not equilibrium metal-loaded folding states or tests of dimerization. The full observed site inventory remains explicit in the embedded path: 8DQ2 has three La and Na; 8FNR four Dy; 8FNS four Nd before replacing target EF3. Exchanging only EF3 therefore creates a conditional mixed-metal state for most endpoints. Spectators are frozen formal monopoles, with no responsive electronic or validated ion-LJ treatment.

## Finite prospective calculation table, no submissions here

1. Finish the current50-atom frozen-f origin and directional derivative gate; distinguish numerical executability from physical accuracy.
2. If passed and a compatible La recipe is fixed, evaluate/reuse **only scientifically identical new-Hamiltonian cells** on all three compact EF3 sources:3sources×2metals=6origins. This is a cheap local-chemistry screen, not the embedded test. Preserve the two Hans-versus-Mex D values separately.
3. For the already prepared complete electronic regions, use3sources×2metals×A/B=12cells, subject to a measured feasibility gate. Report each source's `E_M(B)−E_M(A)` and `delta_env = work_La−work_Dy` independently. Here B moves a donor inside QM, so this is an embedded donor-response test, not the earlier fixed-core external-field perturbation.
4. For both A and B, report `D(8DQ2,8FNS)` and `D(8FNR,8FNS)`. Retain all failures and both source contexts. A sign disagreement is source/state sensitivity; no best-source selection, band refit or unknown-label rescue.

These stages are a design, not automatic execution authorization or a request to repeat the existing scout. Root owns gates, exact Hamiltonian, reuse and manifests. Prior unsuccessful explicit-f and native/GFN2 outputs remain separate evidence.

## Biological evidence boundary

The consumed evidence is recorded in `diagnostics/lanm_series_followup_20260923/REPORT.md` and its superseding `DY_TRANSFER_REPORT.md`: the Hans primary paper (Nature2023, DOI10.1038/s41586-023-05945-5) compares pH5 CD apparent affinities with prior Mex data (Deblonde2020, DOI10.1021/acs.inorgchem.0c01303). The archived interpretation is greater Hans light/heavy discrimination. These are protein-level apparent-affinity comparisons involving folding/cooperativity and Hans dimerization, not a contemporaneous assay batch or independent site-resolved labels. Exact Mex numbers/uncertainties were not recovered in that record; do not invent a quantitative target. No new primary-literature verification or reserved-outcome reading occurred for this note.

The existing compact model's positive Hans8DQ2 EF3 direction reversed on8FNR (archived report), so transfer to both sources is essential, not fresh blind validation. A matching sign on this consumed three-source set would be progress toward a conditional descriptor, not general within-lanthanide affinity validation. No Kd, ΔΔG, occupancy probabilities or protein-level thermodynamic selectivity can be inferred from these electronic exchanges alone.

## Exact reusable artifacts

- Older Hans8DQ2/Mex8FNS: `workspaces/lanm_series_followup_20260923/prepared_v1/endpoints/{Hans_EF3__La,Hans_EF3__Dy,Mex_EF3__La,Mex_EF3__Dy}/core.xyz`.
- Older Hans8FNR: `workspaces/lanm_series_followup_20260923/dy_transfer_v1/endpoints/{Hans_EF3__La,Hans_EF3__Dy}/core.xyz`.
- Complete regions: `workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/{Hans8DQ2,Hans8FNR,Mex8FNS}/INPUTS.json`, with source receipts, exact atom/cap maps, charges, states and water/spectator inventories.
- Preparation checks/reports: `diagnostics/metal_environment_response_20260926/lanm_ef3_preparation/`.

This note used only consumed structures/results and the prepared finite mappings. No new molecular calls, outcome unblinding, file changes to inputs, or Slurm submissions.
