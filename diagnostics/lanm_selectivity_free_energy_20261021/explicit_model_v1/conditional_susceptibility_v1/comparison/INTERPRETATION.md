# Curvature requires separating driving field from response stiffness

The admitted model has **excess added-water susceptibility despite insufficient
zero-field pair dipoles and collective suppression**. A uniform increase in
polarizability is therefore poorly directed for the measured mode. This is a
fixed-geometry, lab-z result against provisional PBE0, not a unique mechanism or
an assessment of every local electric-field direction.

## What was observed

The table uses energy curvature at h=0.001 au for compactness; both steps and both
estimators are retained in the pinned 28-row result. All differences below are
model minus reference, in a0³.

| Quantity | Model | Reference | Difference |
|---|---:|---:|---:|
| La dimer | 18.819716 | 15.472822 | +3.346894 |
| La far trimer | 26.277563 | 22.016524 | +4.261039 |
| Dy dimer | 12.926084 | 13.158107 | −0.232023 |
| Dy far trimer | 20.785375 | 19.880690 | +0.904685 |
| La added-water increment | 7.457848 | 6.543702 | +0.914145 |
| Dy added-water increment | 7.859292 | 6.722583 | +1.136708 |
| Dy−La increment contrast | +0.401444 | +0.178881 | +0.222563 |

Across the four estimates, the La increment excess is 0.912947–0.919817, Dy
1.136291–1.139074, and paired excess 0.218754–0.223344 a0³. The paired sign agrees;
its magnitude is too large relative to this reference. The reference paired
energy/moment estimates span 0.178244–0.183725, so even a passing panel is not
identical beyond printed precision. The numerical gates and rounding bounds
are not physical error bars. In particular, the Dy dimer's −0.232 offset is an
observed deviation, not a calibrated physical failure threshold.

At zero field the existing W9 pair-interaction dipole sums are La 2.916725 versus
3.929194 Debye and Dy 2.729250 versus 3.959697. Their full-shell collective
remainders are −0.809525 versus −1.684749 and −0.711698 versus −1.819874 Debye.
All eight common-basis trimer E3 values are underestimated by 4.617–8.184 kcal/mol;
their collective dipole magnitudes are 58.8–85.6% low, and all four metal orderings
are reversed. These seven-term observables differ from the conditional
trimer-minus-dimer susceptibility measured here: the latter includes the added
water and its other interactions and is not alpha3 or full-shell B.

The new information is curvature, beyond the previously measured first energy
derivative (dipole). It disfavors explaining the zero-field deficits by uniformly
weak electrical susceptibility. It does not rule out a weak or incorrectly
shaped internal driving field, anisotropy, environment-dependent response,
reference error, or a defect in a mode orthogonal to uniform z. CP and common-basis
tests did not remove the energy discrepancy; they do not make PBE0 exact or
uniquely assign its nonadditivity to polarization or charge transfer.

## One concrete candidate: local field plus local inverse susceptibility

Retain native `uind=u`, `uinp=v`, `A=alpha_inverse−T`, and both original exclusion
patterns. Introduce a smooth metal-local field h and symmetric block-diagonal
oxygen stiffness correction H through **one stationary dual functional**:

```
L = K/2 [uᵀ(A+H)v − uᵀ(Ep+h) − vᵀ(Ed+h)]
(A+H)u = Ed+h;     (A+H)v = Ep+h
U* = −K/2 uᵀ(Ep+h)
dU* = K/2 [uᵀ(dA+dH)v − uᵀ(dEp+dh) − vᵀ(dEd+dh)]
```

This functional is stationary, not a joint minimum in independent u and v.
The physical response matrix A+H must be positive definite. For these complete
water fragments Ed=Ep, so at fixed geometry `u=(A+H)^−1(E+h)`. Changing h changes
the equilibrium dipole and energy without changing the linear susceptibility;
changing H changes both. In compatible units the induced curvature is
`Bᵀ(A+H)^−1 B`, where B maps the uniform field onto dipole sites. Thus stronger
internal drive can coexist with lower differential susceptibility. This is an
effective response hypothesis, not proof of a particular microscopic cause.

For each metal M and intact water oxygen i, set r_i to the metal–O distance,
n_i to its unit vector, and freeze the shared radial shape before development:

```
f(r) = [(1+t+t²/3) exp(−t)]² / [(1+5+25/3) exp(−5)]²
t = 2 r/Å                       # original native S2U exponents 4 Å⁻¹
eta_i = sum over other water oxygens j != i of f(r_j)
h_i = a_M f(r_i) Q_M n_i/r_i²    # frozen Q_M=+3; native field units
sigma_i = f(r_i) [b_M + c_M eta_i]
H_i = [exp(sigma_i)−1] n_i n_iᵀ / alpha_O
```

Only metal types 500/501 with water-O type 36 participate. The three amplitudes
per metal, a/b/c, are dimensionless. h has units e/Å² and H has units Å⁻³;
K is the native electrostatic conversion. The shape is normalized at 2.50 Å,
decays smoothly and is not searched against outcomes. No hard coordination
count is used. eta is exactly zero in an isolated metal–water pair, and h/H
are absent without a metal. Intrinsic POLARIZE values, water multipoles,
water–water kernels and the original solvent/protein parameters remain fixed.

The six-amplitude family separates drive, pair stiffness and additional-water
stiffness algebraically. Retain **c_La=c_Dy=0 as the nested four-amplitude model**
within the same prospectively declared development; six amplitudes are not
established as necessary physics. Nonnegative H alone can only reduce the
uniform curvature, so it cannot increase the slightly low Dy dimer value.
Allowing signed stiffness avoids that directional restriction but requires
actual matrix checks; positive local diagonal factors alone are insufficient.

One bounded initial development domain is |a_M|≤0.10, |b_M|,|c_M|≤0.20, with
|sigma_i|≤0.20 on every consumed geometry and derivative displacement. Require
the assembled dimensionless matrix `Abar = alpha^(1/2)(A+H)alpha^(1/2)` to satisfy
`lambda_min > 1e−6 lambda_max > 0`, finite values, and the unchanged native
induction residual ≤1e−8. This is a prospective numerical stability domain,
not a physical acceptance threshold or a claim that a repair lies inside it.
All coefficients start at zero; none is selected here.

## Distinct leverage, required implementation, and development evidence

The pair-stiffness direction overlaps the tested EXCHPOL idea: positive
stiffness alone worsened the zero-field energy contrast. That result must remain
in the evidence. The new candidate adds an independently constrained internal
field and a smooth other-water contribution, rather than repeating EXCHPOL
alone. The D/M probes already separated drive from metal-mediated mutual
coupling, but their fixed ±5% changes were weak, traded near against far E3, and
did not correct metal ordering. H acts on water-site inverse susceptibility;
it is not another scalar metal–O mutual gamma. The proposed radial h may still
be correlated with the old direct-damping response. A consumed Jacobian/rank
assessment must expose that redundancy before treating all parameters as
identifiable; no extrapolation of the old probes establishes a repair.

Required CPU entry points are `polar/kpolar` for the scoped zero-default data,
`induce0a` for both direct fields and the H-corrected initial residual, PCG
operator and preconditioner, and `epolar0e/epolar1e/epolar3a` for identical
stationary energy. Analysis must include the added `−K/2 uᵀh` contribution;
the single-loop path can retain its contraction only when `udirp` includes h
with the original bare-alpha convention. `epolar1a` needs both dh contractions
and the cross-dipole dH term. Differentiate f, eta, the exponential, n_i and
`n_i n_iᵀ`: eta derivatives produce forces on the other water, and angular
derivatives are essential for transverse/rigid-water motion. Preserve permanent
frame torque and virial accumulation, translation/rotation covariance, and
the unchanged external-field permanent-charge/origin convention. For general
dual fields the field derivative uses one half of u+v, not u alone.

The existing consumed constraints can distinguish these directions without
opening reserved data: all 10 energy and 15 force-fixture checks of the original
413-residual development objective; all eight E3 energies and eight mu3 vectors;
the four exact zero-field dimer/far-trimer dipoles; and the four independent
state curvatures here. Keep both field steps/estimators as correlated numerical
diagnostics. The two added-water increments and paired contrast are derived
constraints, not seven independent susceptibility measurements. Retain near/far
and La/Dy outcomes together. A future fit must declare weights, finite bounds
and stopping rules before execution; it must not double-count derived contrasts
or replace the original targets. Sparse trimer terms cannot claim full-shell B.

First numerical qualification would require exact disabled/zero and metal-absent
replay, actual response-matrix stability, energy/analytic-force and torque checks,
and field-energy/moment consistency on these consumed fragments before any
fit. Existing prototype field guards remain unchanged until a separately
reviewed implementation supplies this path. No native patch, build, run, fit,
production change, GPU or periodic release is made by this note.

## Evidence

`INTERPRETATION_INPUTS.json` pins the exact comparison (41dd4eb3), source review,
model/reference results and reviews, zero-field redistribution, trimer result,
and previous sensitivity/EXCHPOL sources. It retains all 28 factual rows.
This interpretation uses saved data and source inspection only, one local
thread and zero new molecular calls. Root retains reference admission and cost
integration.
