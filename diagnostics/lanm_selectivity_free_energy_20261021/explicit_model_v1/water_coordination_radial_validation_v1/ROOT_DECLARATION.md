# Conditional fresh W9 radial validation

Declared21October2026 after submission of joint fit17676, before inspecting any
new fit trajectory or new validation outcome. This release authorizes source and
geometry preparation only. Root owns molecular submission after the stated gate.

## Scientific question and fixed configurations

Does a water candidate that passes the ten consumed joint-development cells also
predict forces and relative work at new nine-water shell radii?

Start from the exact original aquo9_central geometry admitted in17670. Generate
two configurations by translating every intact water rigidly by−0.10Å or+0.10Å
along its own original metal-to-oxygen unit vector. Keep the metal fixed and each
water's internal coordinates/orientation unchanged. Apply identical geometry to
La and Dy: four new cells, all28 atoms. Do not optimize or inspect energies when
preparing them. Preserve deposited source, typing, chemistry and reference recipe.
Record geometry/topology/contact checks and fail any arithmetic/state discrepancy;
no geometric rejection rule based on eventual predicted energies.

These coordinates and reference outcomes must be excluded from17676 and every
fit objective. Register hashes before any scoring. The earlier four orientation
outcomes remain historical regression evidence; they are not fresh validation.

## Fixed comparisons

Compare each new cell to its matching central W9 cell, reusing exact accepted
central model/reference energies and forces. Require all four individual radial
work errors and both Dy-minus-La radial-work errors≤1kcal/mol in magnitude.
Require all six new rigid-force checks (La,Dy,paired for each displacement) to
meet the unchanged category-balanced RMS≤2 and maximum component≤5kcal/mol/Å,
using the existing body helper. Twelve comparisons total. Missing/failed cells
remain unavailable and cannot pass. Report all cells and both signs.

Reference remains the unchanged provisional frozen-f ECP-PBE0-D4 recipe, states,
bases, grids and thresholds used for the existing central W9 references. No
mixed recipe or counterpoise energy substitution. Native analytic forces only;
this transfer experiment does not certify analytic gradients by infinitesimal
energy differences. Preserve native convergence and final-grid warnings.

## Conditional execution and limits

Preparation: zero native/quantum calls, at most one login thread per agent.
Execution may follow only a terminal17676 result passing every old/new development
control and an independent actual selected-ten replay. Freeze its parameters once.
If that gate fails, hold this panel; do not automatically fit its outcomes.

If released later: at most four new fixed-geometry reference invocations on4MPI
per cell, and eight model programs (analyze+analytic testgrad for four cells) in
one four-CPU task. Use standard_expand/defaultmemory/time and exact scientific
sources. No retries, new optimizer, parameter changes or Cartesian finite
differences are included in this declaration. Model/reference arms may run in
parallel after admission. Root must separately bind actual source/config and arm
completion wakes; source preparation is not evidence of execution.

## Interpretation

Passing is limited fresh radial transfer against a provisional reference. It is
not equilibrium hydration, salt/protein transfer, electronic-state qualification
or an exchange-free-energy result. Failure remains retained validation evidence;
any later use for fitting consumes its independence and requires a new declaration.
No automatic combined donor/water promotion follows.
