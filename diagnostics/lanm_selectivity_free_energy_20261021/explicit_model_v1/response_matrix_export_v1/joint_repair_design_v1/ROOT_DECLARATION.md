# Joint local and collective response repair

Root declaration, 22 October 2026 UTC, before candidate implementation or fitting.
The active La/Dy exchange goal authorizes this contained development. Root owns
execution and integration. Existing production and all earlier results stay pinned.

## Question and reason for this model change

Can one conservative interaction model reproduce the consumed water energies,
forces, pair dipoles, full-shell dipoles and electrical response together?

The completed finite screen (RESULT19cb1ea2, independent reviewdfc04b6c,
admission8be74687) tested81 paired settings. All486 response matrices are stable,
but none passes the original coordination-energy check. Its fixed local-field
and metal-stiffness controls have distinct, conflicting effects. La stiffness
improves susceptibility while lowering an already-small pair dipole. The pair
dipole sums need about1D more response, whereas the full W9 interaction dipoles
are much closer to the provisional references. This motivates independent
additional-water screening of the local drive. Dy's dimer and trimer response
offsets also require independent pair and collective stiffness directions.

These observations do not uniquely identify charge transfer, Pauli repulsion or
polarization as the quantum mechanism. The new family is a testable effective
model. The finite grid does not prove failure of the full continuous old family.

## Model scope for implementation

Use two drive and two metal-stiffness amplitudes per metal: eight response
amplitudes total. Distinguish the isolated metal-water field from its change in
the presence of other waters. Distinguish isolated-pair and collective metal
response. The final explicit smooth functions and derivatives must be reviewed
before any candidate molecular calculation.

Preserve the original water model, permanent multipoles, native mutual tensor,
metal polarizability parameters and original damping coefficients. An added
metal response block changes the effective response; it is not a claim that
the original intrinsic polarizabilities are physically exact. Permit a later
joint refit of only the existing metal-O VDWPR radius and epsilon per metal,
four short-range parameters. Do not also refit POLPAIR damping.

Use the stationary dual energy with both native field channels. All forces must
include derivatives of the added fields, stiffness, radial kernels and other-
water environment. Preserve the native local-frame contributions. Moments and
susceptibility must be conjugate to this same energy.

The fixed-geometry diagnostic's unscreened added radial field is singular at
contact and is not suitable for dynamics. The new implementation must regularize
the added field, for example through the fixed native metal-O charge-field
damping factor, and include its derivative. Preserve the old diagnostic as a
distinct result. Require finite short-distance limits, vanishing separated
fragments, rotational/translation consistency and a bounded response multiplier.
These properties do not substitute for actual stability and force checks.

## Frozen data division

DATA_CONTRACT.md (e249cd61) and DATA_CONTRACT.json (3d65804f) bind the exact
consumed inputs. Preserve the10 original energy/force cells and all25 gates.
Use34 literal pair-dipole vectors, four full-shell C vectors and four independent
state susceptibilities for response development. P and B are derived reporting
quantities. The eight E3/mu3 cycles are consumed regression outside this fit,
not fresh validation. The four untouched W9 radial cells and excluded donor/
protein outcomes remain held.

## Prospective objective

Define J = (J_old_EF + J_pair + J_C + J_chi)/4.

- J_old_EF is the unchanged original413-component engineering objective, with
  its four original equally weighted energy/rigid-force blocks and all existing
  scales. Preserve all25 individual gates separately.
- J_pair averages squared Cartesian dipole errors, scaled by1D, within each
  of the four parent contexts, then gives those contexts equal weight.
- J_C averages squared Cartesian errors of the four full-shell C dipoles,
  scaled by1D, equally over parents and components.
- J_chi averages the four state susceptibility squared errors, scaled by1a0^3.
  Fix the reference estimator to the energy curvature at h=0.001au for all four.

These weights and scales are engineering choices, not physical error bars or
inverse-variance weights. Report every block and every individual check. Do not
duplicate P/B, derived susceptibility contrasts, or the other field estimators
inside the objective. Retain all28 field rows descriptively. Do not use E3/mu3
to select coefficients. Existing1kcal/mol energy and2/5kcal/mol/A rigid-force
limits remain unchanged; no physical dipole or susceptibility tolerance is
invented. A small training loss is not model qualification.

## Work released by this declaration

Prepare the explicit mathematical specification, native implementation map,
objective adapter and a finite numerical qualification protocol in disjoint
agent-owned files. Numerical source controls use one local thread each and
saved or synthetic inputs only. No fit or new quantum calculation is released
by source readiness. Root will bind reviewed source, exact parameter domains,
starting point, finite optimization limits and actual numerical qualification
before launching the single declared development fit.

Arendt owns the new response specification and implementation design. Planck
owns independent force/energy integration review. Singer owns the data contract
and objective adapter. Root integrates and submits. The separate already-queued
tighter Dy-acetate reference comparison remains independent.
