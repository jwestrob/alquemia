# Staged periodic native numerical qualification

Frozen prospective contract; no build, native program, factorization or solve is
released by this document. Stage 1 is implemented/reviewed first. Stage 2 needs
its own source review and root release after actual Stage 1 admission. Dependencies
and exact saved authorities are in `DEPENDENCIES.json`.

## State and scope

Use only the prior terminal **e41f** La candidate and its literal 49.324 Å BIO
fixture: 3,980 waters, seven Na, ten Cl, one La, 11,958 positive-alpha sites in
native atom-major order. Preserve composed parameter/key values, intrinsic alpha,
POLPAIR, topology and all mutual scales. The current calibration is independent.
Retain the admitted initialization settings: orthorhombic conducting/TINFOIL,
9 Å real Ewald cutoff, permanent/polarization alpha 0.4 Å⁻¹, FFTW 64³, both spline
orders 5, vdW cutoff/taper 9/8.1 Å without correction, MUTUAL, polar-eps 1e-8,
100 iterations, one native thread. `EWALD-BOUNDARY` remains forbidden because
the native parser interprets it as VACUUM. Restore the two exact composed response
records in the separately reviewed active harness; the operator-only derivative
input from the earlier PME export is not an active-model energy authority.

The compact response kernel remains unchanged through 20 Å, quintic C2 to zero
at 21 Å, and zero outside; the orthogonal shortest side remains at least 44 Å.
All h/H and coordinate/cell derivatives use that same kernel and switch term.
The current saved static certificate applies only to this central candidate,
geometry, alpha, real split and PME state. Its nominal shift is
`tau=1e-6*(U_B+U_R)+delta_odd=0.004903508115932018`.

The numerical questions are native implementation correspondence and finite
energy/gradient consistency. They do not establish physical accuracy, a new
parameterization, stability along an unsampled trajectory or production MD.
`e_app` stays null/unknown. A rigorous whole-application error bound is optional
stronger evidence; it is not a prerequisite to these ordinary-double tests.

## Stage 1: one isolated operator program

Implement a separate opt-in harness around the native common `dfield0/ufield0`
paths and the already qualified reciprocal primitive; do not replace the PCG
solver or construct a dense global operator. Bind the exact central static
certificate and candidate before entering either formerly fatal periodic hook.
The qualification clone may accept only the exact prepared static state; absence,
bad pins or any state change remains fatal. The existing trees/guards are unchanged.

Let `a` be native positive intrinsic polarizabilities, `P=diag(sqrt(a))`, and
`S=I+P H P-P(T_real+T_self)P`, the admitted unshifted symmetric saved matrix.
For normalized inputs `(z,zp)`, set native moment inputs `(Pz,Pzp)` and measure

`A_native(z,zp)=(z-P*ufield0(Pz,Pzp)[1], zp-P*ufield0(Pz,Pzp)[2])`.

The independently assembled prediction is `(S*z,S*zp)+R_packed(z,zp)`. Native
`umutual1` **subtracts** the gathered reciprocal dipfield from the mutual field;
the normalized reciprocal term therefore enters `A_native` with a **plus** sign.
The h correction is added equally to the two direct RHS; H is subtracted once
from each mutual field. Unequal direct fields/solutions remain allowed.

Use the exact ten vectors from pinned `BACKEND.ten_probes(35874,(0,1,2))`:
three uniform Cartesian axes, three metal basis vectors, and four Rademacher
vectors from PCG64 seed 20261022. For each ordered vector z_j, run the three
ordered packed modes `(z_j,0)`, `(0,z_j)`, `(z_j,z_(j+1 mod 10))`.
This is **30 complete mutual-operator applications and 30 separate reciprocal
reference applications**. The complete applications already contain 30 reciprocal
calls: total **60 mutual reciprocal primitive applications**, not 30. Count the
single direct/permanent PME initialization separately. No induced solve, energy,
gradient or optimizer is called in Stage 1.

Snapshot site/alpha/position/order/grid/cell/qfac and B-spline weights immediately
before the applications, and verify equality with the actual initialized state
used throughout. Compare native primitive outputs to an independent sparse
scatter + complex FFT + unchanged stored-qfac + adjoint gather replay. Retain
the native unnormalized backward factor M=64³. Never symmetrize/clamp qfac or
discard packed-channel leakage. Replay native norm/q summaries with the inherited
`1e-12*max(1,abs(replayed))` bound. Require finite, nonnegative stored q, and report
actual conjugate-odd maximum and leakage; expected tiny nonzero odd values do
not become an exact-evenness assertion.

For each packed application independently, set
`b=1e-11*max(1,max_abs_predicted_output)` and require every output component
difference <= b, both for native reciprocal versus FFT replay and for complete
native action versus `S+R_packed`. This inherits the finite PME qualification's
scaled numerical tolerance; it is sampled correspondence, **not** a uniform
norm bound e_app. Record actual errors, scales and both channels. No native
u/v equality, sampled-PCG convergence or parity pass substitutes for the static
stability certificate. Native moments, model state and weights must be unchanged
by the harness. Scratch must be restored/reinitialized before any later use.

One actual `dfield0` initialization may also save before/after hook RHS and h.
Require each channel's added field to match h componentwise under
`128*N*eps*max(1,max_abs_values)` (the inherited arithmetic helper, N=35874).
This is a readback during the same initialization, not another direct-field call.
Save the exact h/beta/environment values used by the action. Source review must
confirm that common hooks execute exactly once per channel/path.

## Stage 2: finite energy/gradient/virial panel, separately gated

One active central program, 32 displaced active programs and three central bypass
programs: **36 energy/gradient native programs**, each with exactly one analytic gradient call
returning full-precision total energy, all 35,874 gradient components and total
virial. There is no internal Cartesian FD or geometry optimization. Three bypasses
are old library/response absent, qualification clone/absent, and clone/zero
amplitudes, on otherwise identical geometry/parameters/runtime. Require complete
energy/gradient/virial tokens to match exactly across the three bypasses.

Use h=1e-4 Å and h/2=5e-5 Å; each of eight directions has q=±h,±h/2. Let L0=1 Å.
All direction definitions use the fixed central state and stay fixed at both steps:

1. Three metal translations: r_M(q)=r_M+q*e_x/e_y/e_z; all other coordinates/cell
   unchanged. Tangent coefficients are the corresponding metal unit vectors.
2. One intact nearest-water rotation. Select minimum nearest-image O–M distance,
   breaking an exact distance tie by ascending native oxygen atom index. Construct
   local water vectors x_O=0 and x_H=minimum_image(r_H-r_O); let c=mean(x_O,x_H1,x_H2),
   u=minimum_image(r_O-r_M), and n=normalize[u×(x_H1-x_H2)]. Fail if the axis is
   undefined/nonfinite; no alternate direction. Define
   r_i(q)=r_i+[Rot(n,q/L0)-I]*(x_i-c) for the three atoms. Thus q=0 preserves the
   literal coordinates even for an originally wrapped molecule. Tangent is
   n×(x_i-c)/L0. Cell and all other coordinates are fixed.
3. One fixed whole-box direction: in native active atom-major coordinates set
   v_j=rademacher_0[j]/sqrt(3), using the same pinned probe definition above;
   r_j(q)=r_j+q*v_j. Cell fixed. This tests a predeclared combination, not every
   Cartesian derivative independently.
4. Three orthogonal diagonal strains: S=e_x e_x^T, e_y e_y^T, e_z e_z^T.
   r_i(q)=r_i+(q/L0)*S*(r_i-r_M), and cell(q)=(I+q*S/L0)*cell.
   The central metal is the fixed anchor. Grid dimensions/orders, Ewald alpha,
   physical cutoffs and coefficients stay fixed, while native weights, real
   geometry, environment h/H and derivatives are rebuilt for each changed state.
   No off-diagonal/nonorthogonal strain is included.

Use the native orthogonal nearest-image convention (ties away from zero).
Preparation must save the exact direction coefficients, atom identities and all
32 immutable coordinate/cell decks before release, without evaluating outcomes.

For coordinate paths, analytic slope is sum_j G_j*v_j, with G=+dE/dr, not force.
For each diagonal strain it is `S:virial/L0` in the native coordinate-gradient
virial convention, including cell/reciprocal terms. FD is `[E(+s)-E(-s)]/(2s)`.
Require each FD/analytic difference <= **0.001 kcal/mol/Å plus literal printing
uncertainty**, and the two FD slopes' difference <= **0.001 plus their summed
printing uncertainties**, inherited from the admitted directional panel.
For an ES token define half-last-decimal-unit uncertainty b. The slope allowance
is `(b_Eplus+b_Eminus)/(2s)+sum_j |v_j|*b_Gj`; for strain replace the last term
by `sum_kl |S_kl|*b_Vkl/L0`. Do not turn this printing allowance into an SCF,
truncation or physical-error estimate. Record native induction residuals and
termination under unchanged polar-eps/iteration limits; no adaptive tightening,
fallback or tolerance adjustment follows a failure.

At the active center additionally retain both direct RHS and both u/v channels,
full-precision polarization energy and correction readbacks. Use the same
stationary dual reduced energy `-electric/(2*dielec)*u·(E_p+h)` with each site's
actual intrinsic alpha and native units. Compare to the native polarization
term with the inherited `1e-8 kcal/mol + arithmetic(N,values)` contraction gate;
do not compare the polarization contraction to the entire bonded+vdW energy.
This readback/replay adds no energy/gradient program. Any extra diagnostic
operator applications must be explicitly counted in the future Stage 2 source
manifest before release; no hidden repeated molecular evaluations are allowed.

All 32 displaced active states require their own nominal state-specific assembly,
PME metadata and sufficient-factor admission before their response solve. The
central job20053 certificate is reusable only on exact state/source correspondence.
Using the existing external-certificate route requires **32 additional native
operator-only initialization programs**, one per displacement, each with one
permanent reciprocal initialization and no induced solve/energy/gradient. Thus
Stage 2 totals **68 native programs (32 initialization +36 energy/gradient)**,
32 saved matrix assemblies and 32 sufficient-factor attempts; central authorities
are reused. The two-stage prospective total is **69 native programs**, including
the one Stage 1 parity program. These are ceilings, not permission to execute.
These independent certificates can run in parallel Slurm allocations under a
separate declared resource budget; a central pass is not propagated along paths.
Any missing, failed or inconclusive certificate leaves the corresponding native
slot untouched. This contract neither releases those factors nor requires a new
rigorous interval calculation. Costs and program/certificate counts stay separate.

## Stop, evidence and later dynamics

Each stage has a fresh immutable output and exact pins for source, build, candidate,
geometry/cell, alpha/map, settings, certificate and root release. Fail on collision,
missing/nonfinite outputs, mismatched pins or numerical gates; preserve raw failure
and all untouched slots. Do not retry or silently alter solver tolerances. Record
link/build costs, native program counts, direct/reciprocal/full-action counts,
gradient calls, certificates, actual threads, scheduler cost and memory separately.

The measured 75.6-second factor is practical for a finite panel but cannot be
assumed economical at every MD step. A later dynamics design needs an explicitly
qualified update/reuse policy with state-change invalidation or a cheaper sufficient
bound; neither is part of this finite contract. No dynamics guard is opened here.
