# Interaction-model checkpoint — 21 October 2026

## Latest result: the La water response needs repair

All48 matched displaced-water calculations complete with exact native replay.
The two displacement sizes agree, and PBE0 energy differences reproduce its
analytic forces within0.021kcal/mol/Å. No central/displaced energy-stage change
was detected. At the existing2.28816Å metal–oxygen distance:

| Metal | AMOEBA radial water force | PBE0-D4, counterpoise | CCSD(T), counterpoise |
| --- | ---: | ---: | ---: |
| La | 12.139 inward | 4.833 outward | 12.463 outward |
| Dy | 32.808 inward | 31.901 inward | 25.833 inward |

Units are kcal/mol/Å. Both quantum methods oppose the model's La motion here;
changing reference method does not rescue that control. Their magnitudes differ,
and their shared frozen-f/basis approximation is not physical truth. Near agreement
of interaction energy at one geometry did not qualify the force response.

The fixed-La/Dy-only fit remains held. Next is a separately declared paired-water
representability test with both metals adjustable and one consistent energy/force
target. Donor reference comparisons17547–17557 run in parallel. No production
model or free-energy prediction is qualified. Actual water-test cost:7,324allocated
CPU-s; no repeated central calculations. Exact result:
`reference_chemistry_v1/water_force_qualification_v1/collection_v1/RESULT.json`
under the workspace.

## What the completed calculations show

All 28 original frozen-f PBE0 configurations, 28 matching MACE/AMOEBA
configurations, and the six central hydrated donor references are available.
The hydrated comparison reuses two aquo8 references and includes 18 new
configurations for each affordable model.

The table gives balanced, fixed-configuration Dy-minus-La electronic-energy
contrasts relative to the specified aqueous cluster. Units are kcal/mol.
Negative values indicate relative Dy stabilization in this cycle, not a measured
binding preference.

| Central hydrated donor | Frozen-f PBE0 | MACE | Transferred AMOEBA |
|---|---:|---:|---:|
| Monodentate acetate | +0.575 | −1.737 | −6.236 |
| Bidentate acetate | −1.974 | −4.566 | −9.812 |
| Acetamide | +1.775 | +1.206 | +9.163 |

MACE is closer to this reference than the transferred AMOEBA model on these
three contrasts. This does not establish either reference accuracy or a
qualified protein-sampling Hamiltonian.

The hydrated fixtures also have longer metal–donor distances than the original
isolated fixtures. Their smaller MACE/reference discrepancy cannot be attributed
to hydration alone. The matched water-removal control is now complete for all six configurations
with each of the three methods. It preserves every surviving coordinate. Neither
comparison is an ensemble or a free energy; the separately reviewed result is
in `matched_water_comparison_v1/collection_v2`.

Independent review reproduced all 21 hydrated contrasts, 42 rigid-ligand energy
derivatives and six force-error summaries. The largest derivative discrepancy
in the arithmetic replay was below 3e−14 kcal/mol/Å. See
[the review](experimental_target_v1/HYDRATED_COLLECTOR_REVIEW.md).

## Reference qualification is progressing

The six-cell canonical CCSD(T) water counterpoise pilot is complete. On its
frozen-f effective-core model, interaction energies are La −89.4674 and
Dy −104.7223 kcal/mol, a paired difference of −15.2549. This is a metal–water pair
interaction comparison, not a water-to-donor exchange or protein affinity.

All declared native output gates passed, including actual ghost charges,
effective-core assignments, correlated electrons and coupled-cluster convergence.
Basis, electronic-state and numerical uncertainties remain separate. The matched five-cell DFT counterpoise calculation is complete, reusing the
two original complexes. PBE0-D4 gives Dy-minus-La -14.7871 kcal/mol, differing
from canonical CCSD(T) by +0.4678; pure electronic PBE0 differs by +0.3399.
Individual PBE0-D4 interactions are 7.90 (La) and 7.43 (Dy) kcal/mol more
attractive, largely canceling in the paired result. This exceeds the declared
0.25 kcal/mol comparison resolution; neither method is established as exact.
See `reference_chemistry_v1/dft_water_cp_v1/final_collection_v1/REPORT.md` in
the workspace. The next declared correlated comparison is donor/water exchange.

## Matched water response now separates the geometry confound

All18 new endpoints are complete. Keeping every surviving metal/donor coordinate
fixed, water addition changes the Dy-minus-La contrast by:

| Donor | ECP-PBE0 | MACE | Transferred AMOEBA |
|---|---:|---:|---:|
| Monodentate acetate | -38.768 | +43.977 | -39.255 |
| Bidentate acetate | -39.094 | -27.892 | -32.659 |
| Acetamide | -45.364 | -26.785 | -46.547 |

Units kcal/mol. Independent raw-energy, state and force arithmetic checks pass.
The MACE environment response disagrees strongly with this reference. Transferred
AMOEBA agrees better on two increments but retains poor hydrated donor/aquo
contrasts and force responses. Neither is thereby qualified. The separate old
versus matched-dry placement comparison uses the same aquo8 cycle throughout.
See [the complete reviewed comparison](matched_water_comparison_v1/REPORT_v2.md).

## Decisions

- The completed28-cell force-target preparation changes the immediate next
  step. The fixed La control fails the declared rigid-body force thresholds:
  water RMS6.163/max16.143 and donor RMS8.776/max39.674kcal/mol/A;
  five of six La radial-work errors exceed1kcal/mol. Dy-only fitting remains
  unsubmitted because it cannot qualify that paired model.
- Independent source/output audit finds no La parameter-transfer, type/class,
  units or native derivative implementation error. Original La interactions
  reproduce the published SI. But its CP-correlated fitting target differs from
  our uncorrected frozen-f PBE0 force targets. The known monoacetate limitation
  remains; the new La discrepancies cannot all be assigned to force-field error.
- The matched two-step water-force comparison is complete as reported above.
  Keep current fit source/data, failed controls and the unopened validation set
  intact; the result supports a fresh paired-water development experiment.
- Preserve both candidate models and the original baseline; promote neither.
- Extend the completed water reference comparison to paired donor/water exchange.
- The saved attribution localizes paired force errors mainly to metal and donor.
  Prepare a bounded parameter repair with new reserved validation configurations;
  shared ligand/water internal distortion must not be fitted into ion terms.
- Keep actual protein donor multipoles in the prepared peptide fragments.
- GPU numerical parity passed job17383: four native polarization tests/182
  assertions and four La/Dy energy/force/finite-difference controls. The largest
  CPU/GPU energy discrepancy is 3.9e-9 kcal/mol and printed force-component
  discrepancy 4.99e-5 kcal/mol/angstrom. Spatial finite differences pass the
  declared1e-3 tolerance. This does not qualify physical parameters or dynamics.
- Numerical dynamics regression17446 passed127 original assertions: three
  trajectory tests and one initialization test. These use upstream test fixtures;
  they do not establish stability of the prepared protein/water Hamiltonian.
- The continuous metal-parameter path17447 passed all20 calls. Four endpoints
  and two CPU/GPU midpoints agree within3.91e-9kcal/mol and4.99e-5kcal/mol/A;
  spatial finite differences pass, and two-step lambda derivatives differ by
  at most5.01e-7kcal/mol. This establishes the finite implementation test only.
- Twelve hydrated radial native references are complete. The next saved
  derivative review and frozen repair pilot test whether the small parameter
  family can reproduce the target response. The reserved24-cell validation
  panel is prepared and its outcomes remain unopened. No candidate is qualified.

## Exact saved results

- `workspaces/lanm_selectivity_free_energy_20261021/chemistry_comparison_v1/ecp_complete_v1/RESULT.json`
- `workspaces/lanm_selectivity_free_energy_20261021/hydrated_comparison_v1/partial_v1/RESULT.json`
  contains the complete 20/20/8 inventory despite its directory name.
- `workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/canonical_water_cp_v2/collection_v1/RESULT.json`

Source pins, native output pins and costs remain in those results and their
task-local execution/submission records. The full goal still requires qualified
interactions, converged ensembles, balanced exchange free energies, assay mapping
and independent experimental validation.
