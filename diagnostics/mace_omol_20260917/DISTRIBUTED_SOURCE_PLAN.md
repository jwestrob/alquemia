# Next source model: one charge/dipole fit with separate spatial validation

Declared after the failed monopole field screen, before fitting or validation
outputs. Active-goal autonomy applies. This is an electric-density representation
experiment toward the same MACE/native-polarization hybrid, not a biological
classifier fit. The failed monopole data remain development evidence.

## Scientific question and exact scope

Can a uniform, compact charge-plus-dipole representation on the existing physical
source atoms reproduce the saved quantum potential and electric field, including
nearby environmental sites, without new SCF or a costly population analysis?

Reuse the8normalized vacuum QM states, CHELPG projection as an initial prior,
actual center-point potentials from explicit_field_short_v1, and native fields
from qm_electric_field_v1. Verify source density, physical IDs, probe order and
state compatibility explicitly; never match by protein name alone. No new
biological cases, geometry, protonation, donors, water, assembly or nuclear force.

Use every recorded physical QM-projection support site, including the metal and
cap anchors; no synthetic cap sites or response variables. Fit one monopole q
and one Cartesian dipole mu per site. These are a distributed approximation
to the whole QM density, not asserted uniquely observable atomic populations.
No quadrupoles, alternate charge schemes or per-protein model selection.

## Fixed fit definition

For source-to-probe displacement r in bohr:

```
phi = sum(q/r + mu dot r/r^3)
E = sum(q*r/r^3 + 3*r*(mu dot r)/r^5 - mu/r^3).
```

Use all actual environment-atom centers as training observations, both the
native potential and all3native field components. Scale field rows by a fixed
length ell=1Å expressed in bohr, giving potential units for every residual.
Use unweighted rows; do not weight proteins, metals or observations by labels.
Represent fitted dipole variables as mu/ell, so they have charge units.

The prior is the old projected CHELPG q with zero additional local dipoles.
Enforce sum(q)=the formal QM endpoint charge through a linear constraint;
this is a new explicitly constrained density fit, not neutralization of a
protein or redistribution of force-field boundary charges. Do not impose an
unverified molecular-dipole origin convention from printed ORCA moments.

In the nullspace of that charge constraint, minimize
`||A*x-y||² + lambda² ||x-x_prior||²`, with
`lambda=1e-4*largest_singular_value(A_nullspace)`.
This fixed regularization suppresses weakly determined coefficients; do not
tune it after either fields or biological ordering. Solve by stable SVD or
an equivalent verified least-squares method, never unregularized normal
equations. Record singular values, effective rank, residuals, coefficient
magnitudes, constraint closure and all actual fitted moments. Large coefficients
remain visible and require scrutiny before coupling to a physical damping model.

## Separate validation and frozen comparisons

For each training position create one observation point translated by exactly
(0.17,−0.11,0.13)Å. No atoms move. This fixed spatial probe set was not used in
the previous field screen or in the fit. Evaluate native potential plus field
there with the same two central-difference spacings0.001/0.0005bohr. Pack the
center and12derivative offsets into one orca_vpot call per endpoint:8newcalls.
No charges may be refit to those outputs. These are spatial representation
checks on consumed proteins, not blind biological validation.

Apply the existing field screen (weighted RMS≤1e−4au OR relative≤0.10) to
all8endpoints and4Ca−La pairs on the separate probes. Retain the1e−6numerical
field and0.01kcal paired-diagonal-resolution checks and1kcal paired U0 flags.
Use the original associated environment-atom polarizabilities as diagnostic
weights only. Also require the historical potential screen RMS≤0.005au OR
relative≤0.10 for endpoints and pairs. Report training and validation separately,
including original/projected monopoles at the identical validation points.
Do not use the old far-only subset as an acceptance route.

Require sum(q) closure≤1e−9e. Verify physical atom order and charge/dipole units,
analytic field signs against differentiation of this representation, and rigid
covariance under the established proper rotation/translation, including vector
rotation of mu. Numerical tolerance for representation replay1e−8atomic units.
Source moments at fixed geometry do not supply valid total nuclear gradients;
mechanical/entropy corrections remain unavailable.

## Execution and stop boundary for this experiment

Eight actual fits and8native potential calls; no DFT/SCF, MACE, FF energies,
new population analysis, optimization or scoring. Use the established shared
8CPU/16GB allocation; one worker/thread per endpoint. Earlier native field
work took203job seconds/711.484CPU seconds; the new probe count is13rather
than12per site. Record fit time, memory, all failures and total execution cost.
No project time/CPU cap. Preserve originals and freeze source/implementation
hashes in workspaces before execution. No fit or validation task has run yet.

A pass qualifies only the fixed-geometry electric representation on these
probes. It does not qualify Tinker source damping, full physical cavity,
covalent boundary charges, matched electrostatic subtraction, gradients or
predictive usefulness. Integrating it into the full model needs those definitions
and its own physical/partition test; no scalar is manufactured here. A failure
is retained without regularization/model/threshold selection to rescue it.
