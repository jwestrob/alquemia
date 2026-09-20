# Same-potential H convergence continuation v1

Frozen2026-09-20 after the three-preparation diagnostic, before this run.
The parent invited this narrow continuation under the overnight authorization.

One additional numerical minimization, no molecular scoring: load the exact
captured System XML and unrounded repaired coordinates from
`workspaces/accommodation_controls_20260920/hydrogen_diagnostic_v1/sample_4/`.
Use the same native OpenMM8.5.1 CPU1 minimizer, tolerance1kJ/mol/nm and a
**500-iteration maximum**. This changes only how far the same minimization
is allowed to proceed. No new H addition, RNG, potential, inventory or state.

Report accepted iterations and final force RMS over **mass-positive movable H
only**, not over fixed-heavy forces. Keep the prior exact-heavy, identical
topology, zero sub-0.45 Å pairs,0.8–1.5 Å H-parent distances and2.5 Å maximum
H displacement guards. Compare to the recorded initial coordinates. Mark
strict mobile-force convergence only if RMS<=1kJ/mol/nm. No adaptive retry.

This is a generic hydrogen-placement potential. Its nonbonded term has no
bond exclusions, so its H-parent bond equilibrium is longer than a physical
C-H/N-H bond. Convergence of this numerical potential is not validation of
DFT hydrogen geometry or an affinity calculation. No scoring/default change.
