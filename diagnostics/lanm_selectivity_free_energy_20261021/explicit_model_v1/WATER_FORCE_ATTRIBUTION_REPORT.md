# Shared internal-force mismatch and metal-dependent response are distinct

The saved analysis completes all eight central cells: three hydrated donor families plus aquo8, each for La and Dy. It reuses the exact AMOEBA, MACE and frozen-f PBE0 arrays. All existing individual and paired component-RMSE values replay, category coverage is complete, and rigid/internal projections pass reconstruction, orthogonality and squared-norm closure. Serial stdlib arithmetic took **0.577 seconds**, with zero molecular calls.

## Large individual error is primarily ligand internal distortion

Fractions below are percentages of the complete configuration's **squared force error**, not energy contributions. Each entry gives La / Dy.

| Hydrated donor | Metal | Ligand | Water |
| --- | ---: | ---: | ---: |
| Monodentate acetate | 3.59 / 1.00 | 80.09 / 83.27 | 16.32 / 15.73 |
| Bidentate acetate | 0.33 / 1.24 | 72.37 / 75.00 | 27.30 / 23.75 |
| Acetamide | 0.10 / 2.51 | 85.56 / 81.31 | 14.34 / 16.18 |

Almost all of each ligand's error remains after its rigid translation and rotation are removed. Its component RMSE is 15.01–24.20 kcal/mol/Å. Thus the individual 8.99–14.10 total RMSE cannot be interpreted as an ion-parameter error alone. The internal projection is a geometric diagnostic: it includes intraligand restoring forces and environment-induced distortion, and does not uniquely separate bonded from nonbonded Hamiltonian terms.

The water component RMSE is 5.27–6.61 kcal/mol/Å, mostly in internal distortion. Even the aquo8 controls have water errors of 5.9167 (La) and 5.9557 (Dy), with internal water distortion contributing 87.26% and 94.12% of the entire aquo configuration's squared error.

## Exact water geometry explains a shared restoring-force contribution

Saved waters have O–H lengths 0.9572 Å and H–O–H angle 104.52°. The executed water parameters have bond minimum 0.9572 Å, angle parameter 108.50°, and Urey-Bradley H–H reference 1.5337 Å; the saved H–H distance is 1.51390065 Å. The angle parameter is a term minimum, not a claim about the isolated or condensed model's total equilibrium angle.

The exact native bond/angle/Urey formulas give water force component RMS values approximately 0, 5.426168 and 0.141869 kcal/mol/Å respectively. Their combined force RMS is **5.369849**. Reconstructed aquo8 bonded scalar components match the saved native component energies within 1e−7 kcal/mol. This was arithmetic using pinned formulas and saved coordinates, not a new native energy evaluation.

That is a substantial shared internal restoring force at the selected water geometry. Similar magnitudes alone do not prove it explains the full model/reference discrepancy: the reference also contains intramolecular and environment-induced forces. No term was deleted, and no water angle or parameter was retuned. At identical La/Dy geometries, these metal-independent water bonded terms cancel exactly in the paired force difference.

## Paired discrepancy remains on the metal and donor

For the AMOEBA error in the matched **Dy-minus-La force difference**:

| Hydrated donor | Metal squared-error fraction | Ligand fraction | Water fraction |
| --- | ---: | ---: | ---: |
| Monodentate acetate | 43.30% | 40.28% | 16.42% |
| Bidentate acetate | 35.87% | 43.93% | 20.19% |
| Acetamide | 43.28% | 53.84% | 2.89% |

Water paired component errors are only 0.809–1.007 kcal/mol/Å. Most of their contribution for acetate is rigid translation; internal-water error contributes just 3.15% / 4.28% of the complete paired squared error. For acetamide, the corresponding internal-water fraction is 1.02%. The large acetamide paired response is therefore not explained by shared water angle/bond terms. Ligand and metal response remains material, with the asymmetric ion/donor parameter provenance described in [HYDRATED_TRANSFER_REVIEW.md](HYDRATED_TRANSFER_REVIEW.md).

## Consequence and limits

The current evidence separates two questions: substantial common ligand/water internal-force differences, and a remaining metal-dependent donor interaction discrepancy. Do not fit ion parameters to the total individual force RMSE or claim that changing the water angle repairs the paired response. The prepared protein fragments permit a subsequent finite qualification with actual protein donor parameters; their scoring remains unrun. No broad MD or automatic parameter fitting follows.

The PBE0 reference retains its electronic-state and derivative-qualification limitations. These are fixed-geometry force comparisons, not physical truth, coordination populations or binding free energies. Projection uses unweighted Cartesian geometry so squared norms match the reported component RMSE; it is not a mass-weighted vibrational analysis. All per-water, per-category, model-force and error decompositions are retained.

## Pins and execution

- [Result](../../../workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/water_force_attribution_v1/RESULT.json): `9f729a13a73b49709a8a1c359c609ea507c1f18a531f681ea8b72d186275c156`.
- [Source](ATTRIBUTE_WATER_FORCES.py): `55681fa8dbe1b3bfbbad125789b1fc81d5549df3ab9397118ea71004dbc2e146`.
- [Declared scope](WATER_FORCE_ATTRIBUTION_SPEC.md); root comparison `e096e732d70a7e352466ab8a1dc46ff19de8427416a0f21e6dcce3a960352c7c`.
- The result pins every parent input and the exact native source formulas. All thread controls were one. Command: `python ATTRIBUTE_WATER_FORCES.py --root workspaces/lanm_selectivity_free_energy_20261021 --output workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/water_force_attribution_v1`, using the source's diagnostic path from the repository root.
