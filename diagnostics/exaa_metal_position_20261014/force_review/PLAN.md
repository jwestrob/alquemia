# Saved forces under the fixed-ion constraint

## Question and finite inputs

Do the archived MACE vacuum forces still exert appreciable force on the fixed
metal after donor accommodation, and does that force point along the tested
reciprocal metal displacement? This is a diagnostic of the proposal Hamiltonian,
not a composite force calculation or an uncertainty estimate.

Use the exact consumed ExaA two-Ca sample0 control and sample3 target from
inventory SHA757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343.
Retain both metals for each original origin, adaptive_Ca and adaptive_La geometry
(12 saved cells), and each transplanted origin/adaptive_La geometry from MACE4130
(8 saved cells). The denominator is20. No new geometry, model evaluation,
optimization, solvent or DFT calculation is authorized by this diagnostic.

## Definitions, fixed before execution

For each host, define u as the unit vector from its original metal position to
its actual admitted transplanted metal position. Use this same u for all five
saved configurations and both metals in that host's Cartesian frame.

- Report the saved metal force vector, its Euclidean norm and F_metal dot u,
  in eV/Å and kcal/mol/Å. Positive projection pushes toward the transplanted
  position; the corresponding local energy derivative is its negative.
- Reuse the archived physical Kinematics/project implementation for all donor
  angular modes, including its cap Jacobian. Report signed energy derivatives
  and physical-heavy-displacement-normalized derivatives, identifying the
  unchanged four modes selected for the old donor search. Do not optimize or
  select new modes using these values.
- Reconstruct each original saved full_q and add only the actual metal
  translation for transplanted cells. Require coordinate replay within the
  existing1e−8 Å preparation tolerance. Verify original own-metal active angular
  projections against their saved gradients within the existing1e−8
  kcal/mol/radian replay tolerance. Caps follow their original mapping; no
  independent cap coordinate or cap force interpretation is introduced.
- Check exact source, state, model, request, coordinate and force-file pins;
  require192 atoms, Ca−3/La−2 and multiplicity1. The physical metal translation
  modes must act only on the true metal, so projected and direct metal forces
  can be compared without an arbitrary cap contribution.
- Missing or incompatible cells remain explicit in the20-row output. A failed
  donor replay does not invalidate an otherwise matched saved metal vector;
  those two diagnostic statuses are separate.

No numerical force threshold will be tuned after inspection. Compare magnitudes,
directions and residual donor loads across the complete saved set. A finite
0.364 Å transplant can cross curvature: endpoint force directions are local
derivatives and cannot by themselves establish an optimal ion position or the
sign of a finite energy difference. Any nonzero residual under a constrained
search is not proof of optimizer failure. Solvent forces, experimental geometry,
thermal populations and class improvement remain untested here.

## Execution

Root submits a single-CPU normal-priority Slurm analysis on the healthy node.
All numerical imports, coordinate/force reading and projections occur inside the
verified allocation with one numerical thread. Freeze the source snapshot and
configuration before execution. Preserve all20 cells and actual job, timing and
failure receipts. No model worker, engine, folding or optimizer is started.

## Potential later reuse, not implementation in this task

The exact archived joint-metal kernel already adds three metal translations to
four frozen donor angles through `adaptive_metal_proposals.MixedProposal` and
`optimize`, with analytic cap-aware gradients and the existing physical bounds.
The proven seven-CPU `CPU_BATCH.WarmCPU` transport returns the same native energy/
force contract. Its authorization layer requires a correctly bound research
manifest with matching settings; the old public `execute()` entrypoint still
expects its historical GPU32CPU allocation and must not be used as a CPU
launcher. A narrow future orchestration adapter could call the existing kernel
with WarmCPU, preserve the chosen origin baseline and common-pool scoring, and
record its own protocol. This diagnostic neither implements nor launches it.
