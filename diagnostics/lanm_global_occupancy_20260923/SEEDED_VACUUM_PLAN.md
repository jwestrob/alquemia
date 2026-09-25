# Numerical recovery of the two failed native vacuum endpoints

2026-09-25. Continue the already approved whole-protein pilot after Jacob's
“two failed” / “proceed”. This is a technical SCF-initialization recovery at the
same two endpoints, not a new geometry, chemical-state or preference experiment.

Job1216564 completed in21141seconds on112 allocated CPUs. Both ALPB endpoints
converged and passed state/charge/strict-tolerance checks. Both vacuum endpoints
failed after500iterations. In the final50iterations, vacuum energy spanned
28.23505276666083Eh for La and0.451218198188144Eh for Dy; final energy changes
were−2.3Eh and+0.000755Eh. These are unusable endpoints. Merely accepting the
last values or raising the iteration limit would not establish convergence.

Exactly two recovery tasks use the converged **same-metal ALPB** native GBW+xtbw
files as electronic initial guesses for their respective vacuum calculations.
Source coordinates, atom order, charge, physical/effective multiplicity, ORCA,
native parameter hashes and state audit match. This extends the existing tested
matching-basename AutoStart mechanism to a different-medium initial guess; actual
`INITIAL GUESS: XTBRESTART` must be printed before recovery is accepted.

The target contains no ALPB keyword. Native GFN2, geometry, charge/spin, temperature,
TolE1e−10, native mixer and MaxIter500 remain fixed. Only NoAutostart is removed
to enable the verified saved initial guess. MaxCore and MPI counts follow the
actual allocated CPU/RAM policy. Two concurrent workers use all CPUs on one
exclusive memory node (56ranks/cell on112CPUs;112/cell on224CPUs), with25% RAM
headroom. No new MACE, DFT, relaxation, physical states or library expansion.

The primary four-cell report is immutable. Recovery has its own manifest/cache
keys and protocol `nikasha_LanM_native_ALPB_seeded_vacuum_recovery_v1`.
Normal completion, strict state/charge audit, unchanged parameters and actual
native restart are required. A single converged initialization is not proof of
a unique electronic ground state. No affinity/classifier claim follows from
these two cells, and rejected MACE relaxations remain rejected.

Fresh preparation and dry-run on the real pinned endpoints passed. An actual
absent-output collection returned two unavailable values, not zero energies.
Production and the original completed solvent outputs remain unchanged.
The batch independently collects even after executor failure; a detached
job-specific reporter writes a report/vault note and emails the actual outcome.
No automatic additional seed, threshold adjustment or molecular resubmission.
