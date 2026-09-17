# Exact short-range MACE evaluation for physical scaffold response

Under the active MACE goal, before new outputs: implement an early exit from
the pinned PolarMACE forward calculation immediately after its trained
`scale_shift` local readout. Sum that readout, exactly as the installed forward
does for `interaction_energy`. Differentiate THIS scalar with respect to the
physical Cartesian coordinates. Do not relabel total MACE forces, invent
charges, include atomic reference energies, or claim a complete Hamiltonian.

The purpose is affordable full-protein mechanical calculations. Direct use of
this component as a discriminator already failed alpha/GGR ordering; that result
is unchanged. No new biological score or threshold is tested here.

Use the pinned medium model, double precision, existing local memory adapter
and finite runner. No installation, training, GPU environment change, numerical
DFT derivative or new DFT. Retain all exact inputs/electronic state records.
This component depends on geometry/species and has no explicit total-charge or
spin response; charge/density output is unavailable, never filled with zero.

Finite validation set: all 14 completed medium whole-protein panel inputs
(five primary La/Ca pairs plus GGR and4MAE rigid rotations), plus all20 archived
GGR physical-displacement inputs used in the completed curvature screen.
Total34 new short-component energy/analytic-force calls. Compare every energy
to its actual saved full-forward interaction_energy: tolerance1e-6eV. Require
the existing rigid gates0.01kcal/mol and0.001eV/Angstrom. For all four GGR
motion/representation blocks, compare short-component analytic projections
with its own centered energy difference at the archived amplitude: tolerance
max(0.02kcal/mol,5% of |h*g|), for La,Ca andR. Negative signs remain.

No repeated model-weight loading optimization in this initial implementation;
keep isolated worker receipts/recovery. Count startup separately from actual
evaluation and record CPU/GPU/memory. Prior small-core batches were about150s
per20calls; full-forward global times were longer because of the blocks this
shortcut excludes. Actual shortcut speed is unknown before measurement.

Next use, conditional on these gates: define J_M(q)=short_full,M(q)−short_core,M(q)
on identical physical atom motions. A conditional mechanical-response model is

    V_M(q) = E_DFT,core,M(0) + g_DFT,core,M^T q
             + 0.5 q^T K_MACE+GB,core,M q + J_M(q)−J_M(0).

Its local gradient retains g_DFT,core + grad(J); its curvature retains cheap
core curvature + curvature(J). Outside coordinates are explicitly fixed, not
eliminated or assigned zero gradients. This measures response from each frozen
baseline geometry; it omits a vertical environmental-energy correction, is not
the rejected total MACE/GB hybrid, and is not a binding free energy. It must be
separately versioned/calibrated if useful. A full response pilot must first pin
source mappings, coupled coordinates, trust region, DFT checks and task counts.
This engine plan itself authorizes no response score or new DFT calculation.
