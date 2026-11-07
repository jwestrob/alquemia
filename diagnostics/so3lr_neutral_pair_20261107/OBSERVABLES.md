# Endpoint observations, declared before dynamics

Use only the three prepared cis/trans systems and their eventual native
trajectories. First run the same analysis on all six actual initial endpoints.
This establishes their starting coordination and hydration, verifies atom
correspondence, and checks that the analysis handles the real system. It is
not an additional potential evaluation or a synthetic test.

For each ligand and metal, retain all ten possible O/N distances, their atom
indices, and the metal distance from their centroid. Count O/N contacts and
water-oxygen contacts at 3.1, 3.3 and 3.5 A. Report distinct acetate arms as
well as oxygen atoms: either oxygen can coordinate, and bidentate acetate is
not two independent arms. A distance count is a diagnostic, not a prescribed
coordination number. Record the first-shell oxygen identities at 3.3 A so
changes can be inspected; boundary crossings alone do not establish exchange.

Retain covalent-bond length extrema for each ligand, tetrahedral-volume signs
at the already verified carbon stereocenters, and distances of both K ions
from the metals and ligand oxygens. Flag ligand bonds outside 0.7–2.0 A and
stereocenter sign changes relative to that start. These are preparation or
chemical-state diagnostics, not fitted selectivity acceptance criteria.

For originally labeled waters, record O–H length extrema and the number of
O–H distances outside 0.7–1.35 A. A proton-transfer candidate requires its
original O–H distance to exceed 1.35 A and a different O/N atom within 1.25 A.
Save the involved atom indices. This is an alert to inspect a reactive event,
not an automatic species assignment or a reason to suppress proton transfer.

For native HDF5, require float64 positions, Cartesian-coordinate metadata,
matching atomic numbers, cell convention, timestep and frame interval. Analyze
only written frames, preserve partial trajectories, and report their actual
duration. Read temperature and timing from native log columns; do not infer
their units from the momentum array. Rounded log energies cannot be used for
cross-state thermodynamic analysis. This analysis makes no model calls.
