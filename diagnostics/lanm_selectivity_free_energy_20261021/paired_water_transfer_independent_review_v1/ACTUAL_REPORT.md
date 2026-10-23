# Actual frozen water-candidate transfer review

**Numerical checks pass; full development transfer fails.** All32 cells completed with48 native programs and1380 internal displaced energies. All10 Cartesian checks pass; the largest component mismatch is 8.4e-07 kcal/mol/Å. Raw native E/forces and the complete comparison replay exactly.

| Group | Baseline failed | Candidate failed |
|---|---:|---:|
| fitted_centers | 4/6 | 0/6 |
| water_radial | 12/15 | 3/15 |
| hydrated_donor | 49/57 | 49/57 |
| dry_baseline | 9/9 | 9/9 |

Water radial energy checks pass, but all three compressed-water force comparisons fail. Hydrated donor transfer remains inadequate. The six dry controls are exactly unchanged: their9 failing rigid checks are pre-existing, not damage from changing water parameters.

These are development results against consumed PBE0-D4/ECP references, not independent validation or physical affinity. No new molecular calls were made for this review.

Machine-readable review: `workspaces/lanm_selectivity_free_energy_20261021/paired_water_transfer_independent_review_v1/ACTUAL_RESULT_REVIEW_v2.json`.
