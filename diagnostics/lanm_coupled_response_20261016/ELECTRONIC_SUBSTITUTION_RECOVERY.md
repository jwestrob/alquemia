# Archived charge-field recovery

Job12770 completed its saved-data traversal but left the full comparison
unavailable: all16 LBT energy-only cross receipts lack the top-level
`charge_closure_e` field used by the gradient receipt schema. The20 other LBT
cells and10 Mex cells passed; the old12736 comparison replay passed. This is a
schema failure, not a failed molecular calculation. Preserve the original
adapter, config, result and partial menus; no full-pool interpretation follows.

The separately versioned V2 reads the actual cross-cell `charges_e`, requires
exact equality to the archived charge file and solver-observation vector, and
recomputes `abs(fsum(charges)-declared_charge)` on Slurm. Both that value and the
absolute recorded signed `solver_observation.charge_sum_residual_e` must be
finite and at most1e-5. The other receipt types retain their original nonnegative
top-level closure check. The exact `--sp` cross-state declaration remains required.

All46 cells,15 menus, source/state/runtime/geometry pins, Hamiltonian, numerical
thresholds and interpretation are unchanged. No new molecular calls or geometry.
V2 binds the exact failed12770 result and requires the original scientific inputs.
Root owns the fresh one-CPU saved-only submission after independent review.
