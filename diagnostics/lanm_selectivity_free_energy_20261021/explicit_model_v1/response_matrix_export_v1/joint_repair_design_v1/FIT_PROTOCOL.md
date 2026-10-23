# One bounded joint development fit

Prospectively declared by root before candidate molecular evaluations or fitting.
This protocol complements ROOT_DECLARATION.md and PARAMETER_DOMAIN.json(a9e53eb5).
It does not release execution before the new native energy/force qualification.

## Inputs and objective

Use the exact12 box parameters in PARAMETER_DOMAIN.json, zero response amplitudes
and selected17676 radius/epsilon values initially. Keep the original native
POLPAIR and all nonpermitted terms fixed. Use objective_v1's declared531-entry
weighted residual and four equal engineering blocks. All25 original energy/
force gates remain separate. There is no assigned physical response tolerance.

Each full vector needs46 unique metal fixtures:10 energy/force cells,34 literal
dimers, and the two W9 far trimers used for the four susceptibility states.
Use analyze EM for46 and analytic testgrad for10, at most56 programs per cold
full vector. Export response matrices only for the four susceptibility aliases.
Reuse the two bare ions,46 metal-free terms and unchanged-metal evaluations
exactly. Reuse the admitted starting values after zero/off numerical parity.
Evaluate the remaining six trimers only after selecting the final candidate;
all eight E3/mu3 cycles remain consumed regression outside the fitting objective.

## Finite algorithm

Normalize every box parameter to[0,1]. Start from the declared vector only.
At each of at most24 outer steps, evaluate clipped normalized±0.005 probes for
all12 coordinates. Record the actual span, including one-sided box boundaries.
Derivative probes are never candidates for selection.

Build the residual Jacobian from those24 probes. Compute an unrestricted
Gauss-Newton direction with NumPy least-squares relative cutoff1e-6 and scale it
to the box and normalized±0.025 trust region. Separately solve the bounded
linear least-squares problem over the same feasible trust box with SciPy BVLS,
tolerance1e-12 and maximum1000 iterations. Record solver status and KKT residual;
an unsuccessful/nonfinite direction is a failure, not permission to substitute
another solver silently. The older3^n exhaustive active-set routine is not used
at12 dimensions. Its small cases can provide independent source controls.

Evaluate four proposals in fixed order: half/full scaled Gauss-Newton, then
half/full bounded least squares. Complete the batch before choosing. A strictly
lower joint objective wins; the current vector wins ties. Persist rejected
trials and exact provenance. Select the final accepted vector by this rule,
not by the reference E3/mu3 outcomes or a favorable metal-preference sign.

Stop on no improvement after a complete four-proposal batch,24 completed outer
steps, a declared resource ceiling, or explicit numerical/native failure.
Passing the25 E/F gates alone does not terminate minimization of the added
response objective. Report the final E/F gates and every response error even if
the joint objective improves. No alternate start or automatic extension follows.

## Resource bounds and reporting

One standard_expand allocation,12 CPUs, scheduler-default memory/time, at most
12 single-thread native workers total. Maximum672 new parameter vectors
(673 including the reused initial vector), and24000 native programs for fitting.
With unchanged-metal reuse, each derivative batch needs at most672 programs
and each four-trial batch224;24 such steps fit below21504 programs. Failed
programs count. Native numerical qualification and postselection regression
have their own finite inventories and costs.

Preserve every molecular input/output and failure; no automatic native retry.
Record total scheduler allocation cost, child CPU time, wall time, cache reuse,
actual vector/program counts, final parameter-bound contacts and last Jacobian
origin. A cap or local stop is not proof of global nonrepresentability.

A completed fit is consumed-data development, not physical qualification or an
exchange free energy. Root reviews the frozen candidate and actual independent
numerical evidence before any untouched validation data are evaluated. The
original validation criteria remain fixed. No donor, protein, periodic or GPU
release is implied.
