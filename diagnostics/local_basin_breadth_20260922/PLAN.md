# Direct one-coordinate basin breadth — frozen pilot v1

## Agreement and question

Jacob approved the three parallel pilots recorded in
`diagnostics/nikasha_parallel_pilots_20260922/PLAN.md`; exact consumed sources are
pinned in its INPUTS.json. This branch asks whether the *width* of a physically
real donor-torsion well contributes information beyond its sampled minimum.
It does not launch the separate PQQ redox proposal or change production.

Four sources: 1H4I, 4MAE, canonical q9z4j7-pqq-la_model, and
`a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4`. Class labels are for
reporting only. Select the first sidechain_torsion in the source's already frozen
common four-mode selector. This gives A/177/chi2, A/301/chi2, A/213/chi1, and
A/290/chi2, respectively. No new force weighting or label-based mode choice.
The final mode is the actual Asn donor torsion, not a carboxylate surrogate.

## Coordinates and measure

Only this one primitive varies; every other full-q entry remains zero. Source
heavy atoms, hydrogens, caps and charges follow the existing exact Kinematics
mapping. Metal, PQQ, waters and spectator coordinates remain fixed. Both metals
share a single physical mode/domain and unchanged paired chemical states.

A rigid group rotates about its fixed source bond axis. Let r_max be the maximum
perpendicular distance of a moving physical heavy atom from that axis. Freeze
`theta = min(0.8, 2 asin(min(1, 0.8/(2 r_max))))` radians. Thus every grid point
obeys the existing 0.8 Å heavy-displacement limit. No energy-dependent domain
shrinkage/extension. Derived domains (radians):

| Source | theta |
|---|---:|
| 1H4I | 0.347769124214565 |
| 4MAE | 0.7492297928256881 |
| Q9Z4J7 | 0.23008816332005527 |
| A0ACD6B9F2 Ca sample4 | 0.7048068578544397 |

Use dq/(2π), with q in radians. In this one fixed-axis rigid rotation, each
physical atom's squared Jacobian norm is constant. Consequently the rotational
mass metric is constant; it is identical across metals because the moving atoms
and their masses are identical and the substituted metal is fixed. Any common
metric factor cancels from the paired free-energy difference. Preparation checks
this atomwise norm invariance and analytic Jacobian on real coordinates. Synthetic
caps belong to the electronic representation, not additional thermal degrees of
freedom. The measure is not generalized to coupled torsions or moving metals.

## Finite calculations, energy and temperature

65 uniform nodes on [-theta,+theta], including q0. The 33-node coarse grid uses
every second fine node. The fixed inner domain [-0.75 theta,+0.75 theta] contains
49 fine/25 coarse nodes. Use composite Simpson weights and stable log summation.
All required nodes remain in the denominator; a missing/invalid node prevents a
qualified correction. Unsupported geometry is never assigned an artificial barrier.

For each node and metal use the existing native float64 MACE OMOL plus native
GFN2 ALPB-minus-vacuum composite:

`U = E_OMOL * EV_TO_KCAL + (E_ALPB - E_vacuum) * HA_TO_KCAL`.

Reuse exact archived q0 cells. Other source pool candidates remain separate and
are not inserted into this quadrature. The shared finite-candidate scoring
adapter provides identical state checking, native runner and MaxIter500 numerical
policy already qualified in the source pools. No orbital restart, DFT, MD,
optimization, surrogate potential, spring or Hessian.

At most 520 metal/geometry cells total; eight exact q0 cells are reusable, leaving
512 new MACE evaluations and 1024 native GFN2 singlepoints before further exact
cache reuse. This is roughly 10^5–2×10^5 allocated CPU-seconds and minutes of warm
GPU work based on prior 240-GFN scoring; it is a scale estimate, not a cap or
measured cost. Record actual finite manifests and receipts, including failures.

Nuclear configurational temperature is fixed at 300 K; RT=R*300/4184 kcal/mol.
The existing GFN electronic smearing is a separate method parameter.

## Quantities and predeclared gates

For each metal, integrate
`Z_rel = integral_D exp[-(U(q)-U(0))/RT] dq/(2π)`;
`F_rel = -RT ln Z_rel`.

Also report `W = integral_D exp[-(U(q)-U_min)/RT] dq` (radians), so
`F_rel = U_min-U(0) - RT ln[W/(2π)]`.

The paired change is Ca minus La, with distinct minima and breadth components:
`deltaR_min = (Umin_Ca-U0_Ca) - (Umin_La-U0_La)`;
`deltaR_width = -RT ln(W_Ca/W_La)`;
`deltaR = Frel_Ca-Frel_La = deltaR_min+deltaR_width`.

Necessary gates, frozen before new energies:

- Complete valid 65-node curves for both metals.
- Each-metal |F65-F33| <=0.05 kcal/mol; paired change <=0.10 kcal/mol.
- Each-metal |F_outer-F_inner| <=0.05 kcal/mol; paired extent change <=0.10.
- Both sampled minima strictly inside the inner 75% domain.
- Each-metal weight outside that inner domain <=1% of the outer integral.

The 0.05 kcal/mol numerical target is ~0.084 RT and is smaller than the existing
0.2 kcal/mol composite numerical scale. Boundary gates prevent reporting a
truncated descending surface as a measured basin. Do not normalize each domain
by its width: that would change the measure during the extent check. Small tail
weight is necessary, not proof of global convergence or of no remote well.

On any gate failure retain curve, raw finite integrals and reasons, but set the
qualified breadth/total correction to null. No post-result threshold/domain
retuning, extra refinement, favorable-point picking or fabricated missing terms.
The outcome is conditional on this one coordinate, frozen spectators and chemical
state. It is not whole-pocket entropy, chemical-state populations, water occupancy
or binding free energy. No inherited classifier threshold or automatic promotion.
