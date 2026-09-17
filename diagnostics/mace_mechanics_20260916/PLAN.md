# Coupled, scaffold-supported MACE/DFT response pilot

Declared before outputs under Jacob's active goal and discretionary pilot
authorization. Prior turn was progress: direct descriptors failed, but MACE+GB
local curvature passed independent archived DFT checks on GGR. The next question
is whether that curvature transfers and can support useful bounded relaxation
when the surrounding protein supplies mechanical constraints.

## Fixed physical systems and two coordinates

Use the existing extended-amide cores for GGR1GLG, alpha1F6S and alpha6IP9.
Also retain the existing connected GGR core for a representation check, never
as a selectable favorable score. Source preparations are the four corresponding
`stage_a_prepared_v1/*/preparation_manifest.json` files in
`workspaces/ggr_mechanism_20260915`. Native r2SCAN-3c/CPCM(Water)/DefGrid3/TightSCF,
singlets and original formal charges; no new chemistry or water selection.

Whole physical systems use the existing five-case global preparation's chainA,
metal and exactly retained site-water inventories/exclusions. This pilot uses
ORIGINAL source hydrogens, reverting the recorded radial H moves in its NEW
workspace only, to match the archived DFT cores. No PQQ case or terminal-OXT
addition occurs in this pilot. Existing H-corrected global data stays immutable.
Record the known stretched-H limitation; no assertion that those hydrogens are
an optimum. Whole protein remains physically uncapped; QM caps have traceable
source links and no independent degrees of freedom.

Choose the nearest selected backbone carbonyl oxygen by source distance, with
source identity as a tie-breaker: GGR Gln140, alpha Asp84 in both structures.
The actual bonded amide neighbors are GGR Ile141 and alpha Leu85. Require the
complete C/O/N/H peptide and both retained C-alpha anchors from connectivity.
Unsupported chemistry or changed donor membership fails explicitly.

Coordinates: (1) metal displacement toward that oxygen, h_s=0.02 Angstrom;
(2) rigid peptide-unit crankshaft about its two source C-alpha anchors,
h_theta=1 degree in radians. Outside physical coordinates are explicitly fixed.
Keep all atom identities/protonation, caps, waters and typed donors unchanged.
Maximum heavy displacement0.05Angstrom; evaluate coupled corners as well as axes.
Initial trust rectangle |s|<=h_s, |theta|<=h_theta. No spring fitting, negative
eigenvalue clipping, entropy or expansion of this rectangle based on a label.

## Model and exact accounting

Medium MACE remains primary. Use the unchanged validated analytic core model
and OBC-II solvent recipe, with a new learned charge distribution at each point.
Use the separately validated exact short-range component/analytic forces for
the whole protein and identical physical core. Define

    C_core,M(q) = MACE_vac,core,M(q) + GB_core,M(q)
    J_M(q) = short_full,M(q) - short_core,M(q)
    V_M(q) = DFT_CPCM,core,M(0) + g_DFT,core,M^T q
             + 0.5 q^T K_Ccore,M q + J_M(q) - J_M(0).

The response gradient is g_DFT,core+grad(J), including the scaffold linear term;
the cheap curvature is K_Ccore+K_J. No assumption of zero exterior forces: exterior
coordinates are constrained, not integrated out. There is no vertical J energy
added to the baseline and no second full solvation energy. This is a conditional
mechanical response descriptor, not the failed direct MACE/GB energy and not a
complete binding free energy. It omits long-range scaffold electrostatics and
does not claim a self-consistent embedded DFT gradient. Short-range J is a
learned component; its adequacy and cap-subtraction sensitivity must be tested.

Use all nine points of the two-coordinate grid z=(s/h_s,theta/h_theta) in
{-1,0,1}^2 for each endpoint. Recover the full2x2 symmetric projected matrix,
including mixed curvature. Repeat the eight noncentral points at half steps to
check cheap numerical convergence; never construct a whole-protein Hessian.
Source/cap mappings carry analytic gradients into these physical coordinates.
Peptide-path curvature includes the curved coordinate, not a Cartesian eigenmode.

## Finite task inventory and reuse

- Core cheap model: four representations x two metals x17 grid points =136
  MACE calls and136GB calls. Reuse20 exact GGR axis/center calls from the completed
  curvature pilot; execute116 newMACE+116newGB calls. No reuse by proteinID alone.
- Whole short component: three physical proteins x two metals x17 points =102
  calls. Add four alpha core-center short-gradient calls; GGR center gradients
  come from the short-engine validation. Total106 newshort-only calls. All core
  short ENERGIES at grid points are already saved with their full MACE calls.
- DFT: reuse20 executed GGR center/axis points. New alpha center analytic
  gradients4; alpha axis displacements16; same-sign coupled corners on all four
  representations16. Total36 newDFT endpoints. No numericalDFT gradients and
  no gradient/Hessian at every cheap point.
- Conditional minimum validation: only stable, converged, physically supported
  models with an unconstrained optimum INSIDE the frozen rectangle are eligible.
  Up to8 additionalDFT single points (four representations x two metals), plus
  matching whole/core short energies at those exact optima. This is validation
  of a prediction fixed before those DFT outputs, not a geometry search. A failed
  gate, outside optimum or changed donor inventory yields unavailable correction.

All cases are consumed development; the two alpha structures are one biological
group. Their qualified affinity direction is distinct from the GGR same-assay
direction and from PQQ functional class. No threshold refitting or claim of
prospective biological validation. Report every case/invalid denominator.

## Frozen gates and decisions

1. Exact paired state/source invariants; short engine exact-energy/force gates
   must pass before response evaluation. DFT analytic components must be present.
   TightSCF center bridge vs archived normal-SCF:0.05kcal/mol per endpoint/R.
2. Local DFT curvature approximation: max(0.005kcal/mol,25% of |DFT even|), for
   La,Ca,R on every axis and same-sign coupled direction. Anchored displacement
   errors <=max(0.02kcal/mol,25% of |DFT even|), as in the prior curvature screen.
   Negative signs remain. No favorable direction selected after scoring.
3. Numerical refinement: maximum |0.5 z^T(A_half-A_full)z| over the unit square
   <=0.01kcal/mol for core, J and combined matrices, each endpoint and R.
   Evaluate the exact extrema of the2x2 quadratic on the square (corners, edge
   stationary points and interior stationary points), not just its diagonals.
4. Stable positive constrained curvature required; no clipping. Preserve soft
   eigenvalues and condition number; require finite solve residual<=1e-10 in
   relative norm. Report trust_region_exceeded instead of a clipped optimum.
5. At eligible optima, compare predicted change against actual
   DFT_core(q*)−DFT_core(0)+J(q*)−J(0), with tolerance
   max(0.02kcal/mol,25% of |DFT_core(q*)−DFT_core(0)−g_DFT,core^Tq*|).
   Keep predictions pinned before validation outputs. Final corrections require
   this validation; then deltaR=deltaCa−deltaLa. No inherited absolute bands.
6. GGR response correction must differ by<=2kcal/mol between extended and
   connected models with the SAME physical protein. This tests response
   representation sensitivity; it does not erase the baseline-core energy gap.

If gates pass, compare unchanged baseline, extended frozen-core baseline, and
response-corrected rawR, separating preparation from response. Report both
alpha−GGR differences under the common extended policy without fitting labels.
If the response does not improve discrimination, retain it as a failed candidate
despite good mechanical tests and continue the larger goal.

## Execution and cost

Use existing runners/isolated environments, finite manifests and receipts.
Short-component evaluation cost is being measured; full gradient core jobs
previously took roughly250s on16ranks.36 new smallDFT tasks at four concurrent
16-rank workers therefore suggest tens of minutes on64CPU, not a routine cost.
Cheap finite grids are a development workload; measure their actual cost before
claiming a feasible routine score. No project time/CPU budget, no long trajectories,
no default change, no shared-job interference. Source inventories and paired input
integrity are preparation; no scorer is enabled merely by writing this plan.
