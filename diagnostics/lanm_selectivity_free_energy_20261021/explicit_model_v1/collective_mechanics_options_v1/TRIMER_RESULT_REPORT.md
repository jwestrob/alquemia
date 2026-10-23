# Two-water coupling exposes the mechanics failure

The collective-response error already appears with one metal and two waters. All eight fixed-geometry tests underestimate the electronic three-body energy, and the model reverses the La/Dy ordering of the collective dipole magnitude in all four matched near/far comparisons. This supplies a concrete mechanics target; it does not identify one unique missing interaction or establish a binding preference.

## Complete fixed-geometry results

E3 is the extra energy beyond all matching one-body and pair terms. The reference is electronic PBE0 in each original full parent basis, without D4. The model is the unchanged consumed AMOEBA-based candidate. Dipoles are common-origin seven-term differences; their magnitudes below are not total molecular dipoles.

| Parent / pair | Metal | DFT E3 | Model E3 | Error | DFT dipole magnitude | Model dipole magnitude |
|---|---|---:|---:|---:|---:|---:|
| W8 / near | La | 12.889 | 8.079 | -4.810 | 0.6824 | 0.2140 |
| W8 / far | La | 11.743 | 5.069 | -6.674 | 0.3058 | 0.0862 |
| W8 / near | Dy | 14.895 | 8.678 | -6.217 | 0.8032 | 0.1919 |
| W8 / far | Dy | 12.310 | 4.597 | -7.713 | 0.3699 | 0.0583 |
| W9 / near | La | 16.806 | 12.190 | -4.617 | 0.8386 | 0.3456 |
| W9 / far | La | 12.384 | 5.079 | -7.305 | 0.3497 | 0.0867 |
| W9 / near | Dy | 19.961 | 13.140 | -6.821 | 0.9780 | 0.3248 |
| W9 / far | Dy | 12.814 | 4.630 | -8.184 | 0.4100 | 0.0590 |

Energies: kcal/mol. Dipoles: Debye. All raw vectors and coefficients are retained in the result JSON.

## What changes between La and Dy

| Parent / pair | DFT Dy minus La E3 | Model Dy minus La E3 |
|---|---:|---:|
| W8 / near | +2.007 | +0.600 |
| W8 / far | +0.568 | -0.471 |
| W9 / near | +3.155 | +0.951 |
| W9 / far | +0.430 | -0.449 |

The energy error is -8.184 to -4.617 kcal/mol. Collective dipole magnitudes are 58.8–85.6% too small. The near-pair energy difference has the correct sign but is too weak; both far-pair energy differences have the wrong sign. DFT gives a larger collective dipole magnitude for Dy in all four contexts, while the model gives a smaller one.

The discrepancy survives both selected water separations (near O–O 2.51/3.00 Å; far 4.68/4.72 Å), with every metal–oxygen distance held at 2.50 Å. Selection was geometry-only and frozen before execution. This is a comparison on a specified surface, not an instruction to change those geometries.

## Consequence for implementation

An energy-only pair correction cannot change this seven-term remainder or its dipoles. The evidence supports developing an electrical-response extension that can independently change the metal driving field and collective coupling. The native engine uses dual induced-dipole solutions; the energy and derivative paths must be changed consistently, with the original water interactions and disabled-mode baseline preserved. A metal-pair direct/mutual separation is a candidate architecture, not a demonstrated fix or a selected fit.

Do not sum these eight sampled E3 values into the full-shell remainder. They establish that the failure is already present at three-body order, not that higher orders are negligible. Quantum E3 is not uniquely polarization or charge transfer. The separate correlated-reference sensitivity branch remains necessary.

## Admission and cost

All 16 new reference and 16 new model cells completed. Native admission 189a30bf verifies exact states, full-parent basis/ghost masks, effective ECP charges, property-file provenance and raw energies/moments. Independent arithmetic f59223b0 reconstructs all eight cycles without repeating the native suite. All 16 reference calculations retain failed printed density/DIIS criteria under the original accepted stopping protocol; numerical admission is not physical-state or reference-accuracy qualification. Printing bounds exclude those uncertainties.

Slurm jobs 18327–18343: 15,284 allocated CPU-seconds (4.25 CPU-hours), with 16 four-rank reference jobs and one four-worker model batch. Root saved collection: 12.108 s wall / 11.263 CPU-s on one login thread, zero new molecular calls. Source and independent review costs remain separately recorded.

Exact results: workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/collective_mechanics_options_v1/collection_v1/actual_v1/RESULT.json. Original inputs, candidate parameters, failed model tests and production remain unchanged.
