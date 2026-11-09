# Chemical integrity audit — 9 November 2026

The completed scout's saved-data analyzer reports a trans600 covalent bond
longer than12 A, while the original structure passed all graph/bond checks.
Do not launch longer trajectories or interpret this as metal selectivity.

First audit the actual21 optimization frames and100 MD frames. Verify atom
order, the first optimization frame against the exact input, the last optimizer
frame against the beginning of dynamics, and checkpoint versus final HDF5
coordinates using the recorded fractional-coordinate convention. Identify the
specific atom pairs, their original chemical bonds, minimum-image distances,
and the first frame in which each crosses the already declared2 A integrity
threshold. Compare both ligands and retain proton-transfer flags separately.

This is a saved-data audit under Slurm. No new model, dynamics, optimizer,
DFT or experimental label is introduced. If a real chemical-state change is
confirmed, permit one bounded static interface check on the exact saved final
coordinates: compare native checkpoint forces with an intact ASE evaluation
at identical charge0/multiplicity6/head5, PBC and numerical settings. Resolve
the source's force-coordinate convention before assigning any discrepancy to
the potential. Keep all failed calculations. Do not add charge constraints,
bond restraints or a different model to rescue the trajectory.
