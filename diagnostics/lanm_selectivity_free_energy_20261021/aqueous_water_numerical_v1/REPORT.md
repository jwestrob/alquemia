# Water parameter transfer into BIO: numerical controls pass

Both fixed aquo8 cells pass. Job17668 completed0:0, four native programs,
zero new finite differences, four allocated CPU-seconds.

| Quantity | La | Dy |
|---|---:|---:|
| Raw BIO−SI energy, kcal/mol | −0.072472010 | −0.072472010 |
| Energy after subtracting known Urey difference, kcal/mol | −1.81e−9 | −1.81e−9 |
| Maximum remaining Cartesian gradient component, kcal/mol/Å | 1.48e−9 | 1.48e−9 |

Every other reported energy component is unchanged to native printed precision.
The Urey component itself agrees with the analytic prediction within8.19e−9kcal/mol.
The declared energy1e−4 and gradient1e−3 limits remain unchanged.

This confirms the finite namespace transfer of the frozen17605 water parameters
while preserving BIO's different water Urey–Bradley equilibrium distance.
Both full aqueous boxes remain prepared-only. Periodic electrostatics, salt
association, flexible-water sampling and protein transfer are not qualified by
this result. The water fit, independent orientations and prior failures remain
unchanged. Independent source and actual raw-output reviews both passed. The actual review
independently reproduces the Urey arithmetic and all 75 gradient components per
metal (104 checks; zero new native calls). This finite mapping branch is complete.

Exact input, output, release, runtime and scheduler records are under
`workspaces/lanm_selectivity_free_energy_20261021/aqueous_water_numerical_v1/run_v1/`.
