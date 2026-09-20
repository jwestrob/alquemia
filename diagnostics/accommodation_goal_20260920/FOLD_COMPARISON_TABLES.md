# Frozen fold-conditioning robustness comparison

Existing 25 PQQ proteins, 250 saved source geometries. All evidence was consumed; structural samples are not independent biological labels.

| Representation | Method | Descriptor | Correct | Wrong | Inconclusive | Unavailable | Denominator |
|---|---|---|---:|---:|---:|---:|---:|
| context | composite | Ca5 | 21 | 0 | 0 | 4 | 25 |
| context | native | Ca5 | 22 | 0 | 0 | 3 | 25 |
| context | composite | La4 | 23 | 0 | 0 | 2 | 25 |
| context | native | La4 | 22 | 1 | 0 | 2 | 25 |
| context | composite | balanced | 20 | 0 | 0 | 5 | 25 |
| context | native | balanced | 21 | 0 | 0 | 4 | 25 |
| context | composite | Ca125 | 106 | 2 | 1 | 16 | 125 |
| context | native | Ca125 | 108 | 2 | 0 | 15 | 125 |
| context | composite | La100_noncanonical | 97 | 0 | 1 | 2 | 100 |
| context | native | La100_noncanonical | 93 | 4 | 1 | 2 | 100 |
| context | composite | all250 | 228 | 2 | 2 | 18 | 250 |
| context | native | all250 | 226 | 6 | 1 | 17 | 250 |
| context | composite | canonical25_replay | 25 | 0 | 0 | 0 | 25 |
| context | native | canonical25_replay | 25 | 0 | 0 | 0 | 25 |
| context | composite | primary225 | 203 | 2 | 2 | 18 | 225 |
| context | native | primary225 | 201 | 6 | 1 | 17 | 225 |
| core | composite | Ca5 | 21 | 0 | 1 | 3 | 25 |
| core | native | Ca5 | 21 | 0 | 1 | 3 | 25 |
| core | composite | La4 | 22 | 0 | 0 | 3 | 25 |
| core | native | La4 | 23 | 0 | 0 | 2 | 25 |
| core | composite | balanced | 20 | 0 | 0 | 5 | 25 |
| core | native | balanced | 21 | 0 | 0 | 4 | 25 |
| core | composite | Ca125 | 100 | 0 | 10 | 15 | 125 |
| core | native | Ca125 | 104 | 0 | 6 | 15 | 125 |
| core | composite | La100_noncanonical | 89 | 0 | 8 | 3 | 100 |
| core | native | La100_noncanonical | 95 | 0 | 3 | 2 | 100 |
| core | composite | all250 | 214 | 0 | 18 | 18 | 250 |
| core | native | all250 | 224 | 0 | 9 | 17 | 250 |
| core | composite | canonical25_replay | 25 | 0 | 0 | 0 | 25 |
| core | native | canonical25_replay | 25 | 0 | 0 | 0 | 25 |
| core | composite | primary225 | 189 | 0 | 18 | 18 | 225 |
| core | native | primary225 | 199 | 0 | 9 | 17 | 225 |

## State accounting

Prepared sources: 233/250.
Preparation failures: {'exact selected PQQ required': 16, 'overlapping expanded atoms/caps': 1}.
Complete protein pools with invariant core state: 21; incomplete pools: 4; complete pools with differing core states: 0.
Proteins with multiple observed context memberships: 21/25.

Bands are frozen from the existing representation-specific calibration. No thresholds were fit. The original 25 canonical source replays are separate from the primary 225 source transfer geometries.

La4/Ca5 are medians of all four noncanonical La-conditioned/all five Ca-conditioned sources per protein. Balanced is the equal mean of those two medians. Missing members make the corresponding descriptor unavailable; no successful-member selection occurs.

Native OMOL, GFN2 vacuum/ALPB contrasts and transfer corrections remain separate in the JSON. Component medians do not add to the median composite in general. Context membership can change with geometry; fixed-core state invariants and context composition changes are explicitly recorded.

This report is an existing-band transfer test, not a newly calibrated ensemble classifier, broad biological validation, thermal population or binding free energy.

Actual collection: `/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/accommodation_goal_20260920/folds_v1/comparison_v1.json`
