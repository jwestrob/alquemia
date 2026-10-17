# Compare admitted archived partition forces

This finite saved-result analysis uses the exact four previously admitted POINT
records from recovered pool13013 as its force authority: 7CCO/cartesian_01,
baseline/expanded × La/Dy. Compare them to the before-La/before-Dy whole-peptide
reference13424 on the same coordinates, with the original frozen two projections.

This is distinct from failed raw-reconstruction checks13446/13455/13460. Preserve
those failures. Diagnostic13460 finds exact stored-component/full-gradient sums
and exact MM/ALPB replay; cap-scatter MACE/GFN components differ by at most
3.552713678800501e-15 kcal/mol/Angstrom. The precise runtime cause is unproven.
No tolerance is added and the raw bitwise check remains failed. This comparison
does not independently qualify the original analytic gradient.

Require all original source, physical-ID/coordinate, state, common raw ALPB and
native-energy lineage gates. Require finite256 archived component gradients,
exact sum to the archived full gradient, and exact replay of the original baseline
force tables and projections. Bind the four raw diagnostic records to their exact
POINTs. Expanded-after remains unavailable. Retain all four cases, full per-atom
forces and components, and both original frozen projections.

One Slurm CPU, unchanged NumPy2.4.6 environment and thread limits. No new molecular
calls, optimization, geometry, model parameters, fitted thresholds or production
change. Interpretation is same-geometry partition sensitivity; no unique MM fault,
expanded-path validation or preference prediction. Root owns execution and report.
