# Independent review of the completed proposal experiment

Read-only review by the second-shell/nonlinear agent,2026-09-20. Source:
`workspaces/accommodation_nonlinear_20260920/proposals_v1/final_1204171.json`,
its pinned manifest, all60 native proposal receipts and all120 GFN2 tasks.
No new molecular calculation or change to the selection rule.

All60 proposal coordinates exactly reproduce the physical Kinematics mapping at
the recorded q (maximum coordinate discrepancy0Å). Inactive coordinates remain
zero, all angles respect±0.8rad, and source element inventories are unchanged.
The actual low-level inputs exactly match the declared primary native GFN2
recipe at the manifest charge/spin; both media use the proposal coordinates.
Every charge sanity audit passes. Fresh native MACE q0 replay is exactly0kcal/mol
different from its archived value.

Independent direct energy algebra reproduces every endpoint work and every
q0/proposal/selected Ca−La contrast exactly. All49 proposal selections and11
origin selections obey the frozen requirement that composite energy decrease
exceeds0.10kcal/mol. Retaining q0 is an explicit two-candidate energy selection;
no failed solvent endpoint is substituted. All120 proposal GFN2 results exist.

Counts reproduce the reported outcome:25canonical origins are correct under
the existing bands; selected geometries give24correct and1inconclusive(Q9Z4J7),
with no opposite-class call. Consumed crystals retain3/3. The canonical raw gap
increases5.076366→5.844397model kcal/mol, while the old-band Q9Z4J7 call changes
because the numerical interval moves. These facts do not justify threshold
padding or treating improved raw separation as a fully validated new classifier.

PLM labels remain unknown.8344 stays Ca-like;07ab moves to inconclusive under
developmental old-band transfer. Comparison against the full-composite optimizer
must retain the different search objective and primary-versus-Tight numerical
policy; its eight final candidates remain unqualified as complete pairs. No
unconstrained minimum, physical population, DFT agreement or broad accuracy is
established by the proposal selection itself.

One initial reviewer scratch assertion demanded charge sums within1e−8e,
stricter than the existing charge audit; the largest observed discrepancy was
2.3753e−8e. The correct review uses the declared charge-sanity policy. This was
a review assertion error, not a solver or scientific-data failure, and no input,
output, threshold or existing check was changed.

**Conclusion:** no mapping, state, sign, unit, selection or denominator issue was
found that changes the reported conclusions. The potential utility is cheaper
physical proposal generation; established-label fold robustness remains the
necessary accuracy test before any promotion.
