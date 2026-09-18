# Full three-dimensional metal response — declared 2026-09-18

Question: did the previous one-dimensional metal coordinate miss a useful
structural response? In archived analytic DFT gradients, that coordinate
captures only 1.5–31% of the gradient norm for the two alpha-lactalbumin
structures. This motivates testing the complete physical translation subspace,
not a favorable direction selected from a biological label.

Under the active autonomous discriminator goal, use the same four already
consumed preparations: GGR extended and connected, alpha1F6S and alpha6IP9.
Keep source heavy atoms, original hydrogens, protonation, waters, core membership,
formal charges, native r2SCAN-3c/CPCM and MACE-POLAR-medium/OBC2 unchanged.
Original-source stretched hydrogens remain a limitation. The protein and all
caps are fixed. The sole physical motion is the selected metal in Cartesian
x/y/z; no fictitious cap or disconnected-fragment modes, no entropy.

Reuse the previous mechanics center DFT energies/analytic gradients, MACE core,
GB core and full/core short-component results. At metal displacements u, define

    C(u) = MACE_vac,core(u) + GB_core(u)
    J(u) = short_full(u) - short_core(u)
    V(u)-V(0) = [g_DFT(0)+g_J(0)]^T u + 1/2 u^T [K_C+K_J] u.

This is the previous conditional mechanics model in a different, complete
metal-translation subspace. It is not the gradient of the failed conductor
hybrid, self-consistent QM/MM, an absolute affinity or a complete free energy.
Do not combine the dielectric and mechanical corrections.

## Inventory and frozen criteria

Use center, six signed axes, and twelve signed two-axis corners at h=0.02 A;
repeat all 18 noncentral points at h=0.01 A. Thus 37 points/state. Reuse all
8 core and 6 whole centers. New calls:288 core MACE,288 core GB,216 whole short
component. The two GGR representations share the exact same whole protein.
Every paired geometry and typed donor inventory must remain identical.

Construct full symmetric 3x3 finite-difference Hessians; analytic short forces
supply g_J. Cheap finite differences are allowed; no numerical DFT gradients.
Require positive eigenvalues at both step sizes, relative solve residual<=1e-10.
Do not clip negative eigenvalues or remove directions. Refinement error in
predicted energy at the proposed point must be<=0.02 kcal/mol for each C,J,total,
and displacement change between step-size predictions<=0.01 A.

The new protocol's maximum metal displacement is0.20 A (about8% of a typical
2.5 A metal–oxygen distance). This exceeds the finite-difference step but remains
a small local coordination change. It is a NEW validation domain, not a
retroactive change to the previous two-coordinate pilot. All donors must remain
unchanged and direct DFT validation is mandatory; a local Hessian alone does not
justify extrapolation. Outside optima and unstable matrices are unavailable,
not clipped into an apparently valid relaxation correction.

At every stable inside-domain optimum, prepare one native analytic DFT gradient
and matching core/full short evaluations (at most8 DFT and16 short calls).
Freeze coordinates/predictions before outputs. Require actual energy-change
error<=max(0.05 kcal/mol,25% of the nonlinear DFT energy change), lower actual
DFT+J energy, and residual gradient norm<=max(2 kcal/mol/A,25% of original norm).
Report actual energy changes and usefulness alongside failed physical gates;
only paired passing endpoints qualify for a response correction. GGR correction
partition sensitivity must be<=2 kcal/mol. No absolute reference/bands/refit.

Primary usefulness: both alpha structures should have higher corrected R=Ca-La
than GGR under the shared extended-core protocol. Report connected GGR as a
representation check, not a selectable result. These are two biological groups,
not four independent labels. No new PQQ inference is authorized by this plan.

Use existing pinned workers/environments and independent manifests, preserving
all previous failed results and concurrent jobs. This development grid has
hundreds of cheap calls; measure actual receipts before claiming production
cost. No project compute/time cap; no production promotion or baseline change.

Numerical convention frozen before calls: use the half-step matrix for every
prediction. After passing validation, use the actual second-endpoint DFT+J
energy change for the primary response descriptor, retaining the quadratic
prediction separately. Do not select the more favorable one.
