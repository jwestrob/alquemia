# Preparation notes

All four sources prepared successfully. Final immutable run is `prepared_v3`.
The first copy stopped before any energy because the older 4MAE DFT receipt
records a relocated runner path with the same implementation hash. The existing
`hydration_square.endpoint` validator verifies its actual output, original XYZ,
input, SCF convergence, ORCA version and runtime artifacts; those two unchanged
DFT centers are reused. No valid chemistry was rerun.

`prepared_v2` was completed but unsubmitted. Version 3 incorporates explicit
native OMOL head/unit and scalar-only output checks before execution. The
physical states, profiles and finite task identities are identical.

Context sizes are 154 (1H4I), 202 (4MAE), 190 and 168 (the two PLM cases).
Maximum physical heavy displacement is 0.441492 Å. All 64 geometries pass exact
paired-state, source-bond/Jacobian and severe-overlap checks. Four real-fixture
tests pass in 4.666 seconds.

The PLM Asp positive rotation relieves the pre-existing compression without
changing atom inventory: the first site's nearest oxygen moves from 1.766 Å
to 1.925 Å at +0.2 and 2.112 Å at +0.4; the second from 1.847 to 1.998 and
2.183 Å. Opposite rotations compress further and remain in the denominator.
These are fixed geometric interventions, not claimed relaxed structures.

The warm scalar implementation reproduces all four archived control OMOL
energies exactly at saved precision. This verifies scalar reuse but does not
validate new PLM labels or the finite-displacement energy surface.
