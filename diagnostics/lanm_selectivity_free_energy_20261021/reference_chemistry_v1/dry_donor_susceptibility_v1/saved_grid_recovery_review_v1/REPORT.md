# Saved grid recovery evidence

The installed ORCA workflow has executed the required grid controls successfully, but **no successful numerical repair is established**. The three old Mex calculations admitted their native results and printed identical late-SCF and final COSX grids. Their completed derivative comparison still failed. This is evidence for the exact input syntax and its activation, not validation of a replacement energy/property surface.

## Exact input change

For all ten new dry-donor cells, add only:

```text
%method
 AngularGridX 2,4,4
 IntAccX -1,-1,-1
end
```

The actual Mex inputs also contain `DoEQ false`, which predates that grid change. The current donor inputs have no `%method` block or `DoEQ` setting. Do not copy `DoEQ false`. Retain the donor's original `MaxIter 300`; the old Mex admission helper hardcodes 125 and cannot be reused unchanged.

The intended change makes COSX grid2 use the existing grid3 convention: angular level4/Lebedev302 and IntAcc4.871, replacing level3/Lebedev194 and IntAcc4.338. Grid1, final grid3 and the XC grid remain fixed. Final grid3 reintegration remains enabled; matching grids does not mean deleting or subtracting that stage.

## Raw evidence and prospective native checks

[EVIDENCE.json](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_donor_susceptibility_v1/saved_grid_recovery_review_v1/EVIDENCE.json) is available at the absolute workspace path recorded in READY.json. It pins the extraction source, all three old Mex inputs/outputs, the failed historical collector, and all ten original donor inputs/outputs. Each donor record contains its complete expected new physical grid fingerprint. Runtime timing and batch partitioning are excluded from that fingerprint.

| Metal | Unchanged XC points | Unchanged grid1 points | Old grid2 points | Required new grid2 = unchanged grid3 points |
|---|---:|---:|---:|---:|
| La |107266|16780|28597|52712|
| Dy |107343|16804|28624|52765|

Verify the actual angular/radial grids, IntAcc, pruning, weights, cutoffs, partial contraction, H/He reduction, rotational-invariance setting, fitting setting and point counts against those fingerprints. Identical point counts alone are insufficient. All five original fields for each metal have the same physical fingerprints, so the expectations require no field-dependent choice.

In the old Mex central/plus/minus outputs, grid2 and grid3 each contain 1,026,059 points and match in every retained physical setting. Nevertheless, the embedded PBE0 directional discrepancy is −0.123371893 kcal/mol/Å against an allowed magnitude of 0.005379702. Every cell retains failed printed SCF residual criteria; the historical report also retains a minus-cell scalar COSX closure issue. The grid change reduced the prior discrepancy but did not qualify that branch. The historical central additionally differed in MPI count, so its improvement cannot be attributed solely to the grids. See the [preserved historical report](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_coupled_mechanics_20261017/mex_dft/consistent_grid_v1/final_collection_v1/REPORT.md).

## Minimum prospective recovery

Use one uniform new ten-cell numerical recipe: **Dy zero and four fields first; La zero and four fields only after actual Dy native, six-check and cost review**. Keep the original PBE0-D4, basis/ECP, geometries, charge/state, field direction and ±0.001/±0.0005 au values, origins, SCF stopping settings, four MPI ranks, resource rules, property parsing and all sixteen numerical tolerances. No old/new recipe mixture enters the paired result.

The original ten cells remain numerically unqualified: ten native admissions, fourteen of sixteen response checks, both Dy odd-energy/zero-dipole checks failed. Singer's saved diagnosis attributes about 74%/92% of those two discrepancies to final integration, descriptively. Neither that decomposition nor the old Mex syntax evidence permits subtracting a correction, averaging recipes, changing a tolerance, or promising that this recovery will pass. Retain both energy estimators, raw SCF flags, all property/origin/D4 joins and all scheduler costs.

This review used saved files only: zero molecular calls, builds, submissions or producer edits. Root owns any execution. The numerical test does not establish physical accuracy or remove the frozen-f and SCF-state limitations.
