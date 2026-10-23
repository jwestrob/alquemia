# Dry mono-acetate source and preparation review

The declared 16-SP diagnostic is coherent. Independent saved inspection passes
153 input, geometry, state and raw-center checks. The two actual `.engrad`
files reproduce every consumed central force exactly with
`F=-gradient*627.509474/0.529177210903`. The ligand radial derivative is
`-sum(F)·u`; the fixed-axis centroid rotation gives `-tau·n/(1 Å)`, matching
right-handed Rodrigues displacement with `theta=q/(1 Å)`.

| Central derivative (kcal/mol/Å) | La | Dy |
|---|---:|---:|
| Radial | 99.580176972 | 105.901221342 |
| Rotation | 28.341183435 | 26.604998856 |

All 16 literal geometries, exact original decks with only EnGrad removed,
40-electron/189-function inventories, fixed metal and rigid ligand distances
join the preparation. The remaining 12 executions require the four admitted
La-radial outputs and all three unchanged 0.1-resolution scout checks. Full
collection preserves all 18 primary checks and missing/failed slots; numerical
failure cannot become a release merely because native programs terminate.

Both centers print post-SCF gridx3 exchange corrections, D4 and Cartesian
gradients. Raw final energy equals `.engrad` energy; post-integration electronic
energy plus D4 closes. Original Total+1el-Energy stopping, residual flags,
COSX stages, DFT/XC-grid settings and PostSCFGGA are retained. These observations
support testing the saved forces, **not an advance claim that EnGrad and SP
are the same numerical surface**. Actual two-step finite differences decide
only the declared sampled consistency. No center rerun, method change, tolerance
relaxation, molecular execution or held-out result read occurred in review.

Final source/preparation/control/readiness bindings are in REVIEW.json.
Author controls are read, not rerun. A reviewer-only initial numeric-token regex
incorrectly matched gradient timing ellipses as an energy; its source is
preserved as `REVIEW_PREPARED_before_token_fix.py`. The corrected independent
reader requires a numeric token and passes; production source was unaffected.
