# Coupled donor/metal response: fixed path pilot

Declared 2026-09-18 before new energy output, within Jacob's approved autonomous
MACE discriminator goal. Preserve all old results and the production default.
Question: does allowing all declared donors to respond with the metal improve
structural robustness beyond the failed metal-only correction? Use all eight
representations / sixteen endpoints from prepared_v1/preparation.json: both GGR
cores on 1GLG/2FW0/2FVY, and alpha1F6S/6IP9. Two consumed biological groups.

## Model and selection fixed before calculation

Use the same normalized H, vacuum native r2SCAN-3c and masked OMOL descriptor.
Same source microstates, assembly, waters and caps; water/exterior remain fixed.
Physical modes are exactly those in PLAN.md, with all chi bonds and all peptide
crankshafts. No donor selection from favorable scores.

For each endpoint let D be the native core DFT energy, T the masked learned
scalar, x(q) the source-defined physical kinematics, and d=gD(0)-gTcore(0).
The approximate change along the physical path is

    F(q)-F(0) = Tfull(xfull(q))-Tfull(xfull(0))
                + d . [xcore(q)-xcore(0)].

At the center this has the exact saved hybrid gradient. It replaces the
unknown nonlinear DFT-minus-learned core response by its Cartesian tangent;
nonlinear learned full response and exact curved-coordinate mapping remain.
This is neither an exact QM/MM potential nor a solvation/free-energy model.
No CPCM/GB/entropy is added. Curvature comes from the learned scalar, not fitted
springs. Actual native points will test the tangent approximation.

Define the initial downhill direction v=-solve(Jheavy^T Jheavy,gH), using the
unweighted Euclidean physical-heavy-atom displacement measure. This is a
kinematic descent metric, NOT a stiffness or covariance. Require positive,
well-conditioned Gram matrix (smallest/largest eigenvalue >1e-10); do not drop
modes or regularize it to obtain a result. Normalize v so its maximum linearized
heavy displacement is 0.20 A, with every angle <=0.20 radian. If exact kinematics
exceed 0.20 A at any prescribed path point, uniformly shrink this geometrical
amplitude by bisection before any energy evaluation. This depends only on
geometry and the already saved gradient, never on affinity labels or new scores.

Evaluate only lambda=0.25,0.5,0.75,1 along q=lambda*v. Center lambda=0 is reused.
Choose the lowest F among these FIVE fixed points, breaking ties toward smaller
lambda. This is a bounded, discrete response descriptor, NOT a stationary minimum
or converged relaxation estimate. No path extension, fitted threshold or second
round of energy-dependent line refinement belongs to this version.

Reapply the established full-source typed donor selector on every point.
If any changes membership, flag that endpoint unsupported and do not score it;
retain it in denominators. Reject bond changes >1e-9 A or fixed exterior/water
changes. Validate analytic geometry Jacobians at the actual path points with
pure geometry finite differences <=1e-7 A/unit. Caps follow both real anchors.

## Actual validation and prediction criteria

At each nonzero selected point run ONE native core analytic DFT endpoint and
matching learned full/core analytic gradients. Reuse the center if selected.
No numerical DFT gradients or optimization. Actual response is

    deltaH = D(q)-D(0) + Tfull(q)-Tfull(0) - Tcore(q)+Tcore(0).

Retain separate native/core/full contributions and raw Ca-minus-La response.
Native energy-prediction error must be <=max(0.05,0.25*abs(D(q)-D(0)-gD(0).dxcore))
kcal-equivalent, retaining the prior native energy criterion. Actual deltaH must
be negative for a proposed nonzero descent. Compare predicted/actual gradients
on ALL physical modes at q, using each center-defined 0.02 A probe scale:
absolute error <=max(0.04,0.25*abs(actual projected probe)). The 0.04 absolute
floor corresponds to the existing 2 kcal/A gradient floor over a 0.02 A probe;
this tests tangent reliability, not stationarity. Full learned energy replay
must agree with the selected grid energy within 0.01 kcal-equivalent.

Require both response and final-score GGR partition differences <=2 kcal-equivalent
within each structure. All twelve alpha-minus-GGR raw margins must exceed 0.02
for structural-robustness success; missing/unsupported cases remain failures in
that denominator. Report partial changes and all failures. No bands, universal
zero, biological independence of structural replicas, entropy or affinity claim.
All qualified corrections stay null if their endpoint/partition gates fail.

## Finite execution and cost

At most 64 learned full-system scalar calls, then 32 learned gradients and 16
native core DFT gradients at selected points. Previously measured 72-scalar GGR
batch ~1513 GPU-allocation seconds; 8 native DFT gradients ~873 wall seconds on
64 CPU, with four 16-rank slots. Those are engineering scale estimates, not
scheduler estimates or measured current cost. Record actual costs/failures.
No project CPU/time budget. Use existing A5000/16CPU/64474MiB and native quantum
runners, scientific cache keys, locks and source pins. No backend change.
Prepare all source/direction/grid definitions before launching. Use one initial
manifest per physical protein (two GGR representations where applicable), then
one native manifest per engine. No new GPU architecture or persistent-worker
engineering is required for this question.

Deliver actual paired results, physical gates, all denominators, cost receipts,
code/tests, current status and a vault note. A failure rejects this path candidate;
it does not authorize changing its scientific definition after inspection.
