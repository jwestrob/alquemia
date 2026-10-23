# Dry acetate mechanics: matched reference and model derivatives

Both methods pass all18 internal energy/force consistency checks. Nevertheless, the frozen donor model gives the opposite centroid-rotation energy slope from PBE0-D4 for both metals. This is a resolved model–reference disagreement in these sampled directions.

## Matched slopes

Units kcal/mol/Å; rotation uses q=1Å×theta. Positive means energy rises with the declared positive displacement. The smaller0.005 step was prescribed in advance and both-step checks pass.

| Motion | Comparison | PBE0-D4 | Frozen model | Model − reference |
|---|---|---:|---:|---:|
| radial | La | 99.557065 | 98.183174 | -1.373891 |
| radial | Dy | 105.918951 | 119.513365 | 13.594414 |
| radial | paired | 6.361886 | 21.330191 | 14.968305 |
| rotation | La | 28.357940 | -20.390410 | -48.748350 |
| rotation | Dy | 26.601243 | -7.316093 | -33.917336 |
| rotation | paired | -1.756697 | 13.074317 | 14.831014 |

## What supplies the model difference

| Motion | Dy−La van der Waals | Dy−La permanent multipoles | Dy−La polarization | Total |
|---|---:|---:|---:|---:|
| radial | 14.709231 | 0.000000 | 6.620960 | 21.330191 |
| rotation | 7.147166 | 0.000000 | 5.927152 | 13.074317 |

Intraligand derivatives vanish at printed precision for these rigid motions. Model component sums close within their original printed precision. The metal-dependent slopes come from van der Waals and polarization; this partition does not assign a unique quantum cause or prove that either term alone should be changed. Centroid rotation also changes individual metal–donor distances, so it is not a pure angular probe.

The original donor fit had already stopped without a useful local improvement in its existing eight-parameter family. These results motivate testing the donor electrical response alongside short-range mechanics; they do not justify repeating that unchanged fit. The ongoing water-only correction has exact dry bypass and cannot itself repair these dry-force discrepancies.

## Scope and cost

Reference:16 new SPs, two old central gradients reused; independent review6f079b8e. Maximum derivative discrepancy0.040842 and step difference0.032964, original resolution0.1. Energy Check stopping, failed printed density/DIIS/orbital-rotation flags, and grid caveats remain; physical state and reference-method accuracy are not certified.

Model:16 analyze programs, two old central model forces reused; independent review4a0d86de. Maximum derivative discrepancy0.006716 and step difference0.005038, same resolution0.1. Exact17621 selected parameters/binary are retained.

Scheduler cost: reference9916 plus collectors5 CPU-seconds; model11 CPU-seconds. No center reruns, fresh held-out data, parameter change or production promotion. This is a local mechanics diagnostic, not an exchange free energy.

Exact results, reviews, formulas and both-step components are bound in [comparison INPUTS](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/dry_acetate_component_derivatives_v1/comparison_v1/INPUTS.json) and [comparison RESULT](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_donor_fit_v1/dry_acetate_component_derivatives_v1/comparison_v1/RESULT.json).
