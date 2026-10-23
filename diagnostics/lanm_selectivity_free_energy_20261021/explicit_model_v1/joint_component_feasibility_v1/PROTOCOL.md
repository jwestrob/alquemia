# Can the existing water-response model improve the components together?

Root declaration, 22 October 2026, before execution of this diagnostic.

## Question and scope

The completed joint water fit improved a shell-size energy difference through
compensation between pair and collective errors. Does its existing 12-parameter
family have a local direction that improves both components while retaining its
other fitted response properties? This finite calculation examines the saved
response probes and trials. It does not perform a new fit or molecular evaluation.

The active La/Dy goal authorizes this method-development diagnostic. The older
infrastructure-blocked status is superseded. The scientific hold on fresh
validation and production remains: failure of consumed regression is a reason
to investigate, not permission to evaluate untouched validation outcomes.

## Fixed evidence and definitions

Use all 673 complete vectors, all 24 original derivative origins, their 24 probes
per origin and the four pre-existing trials per origin from joint_repair_design_v1/
fit_v1/run_v2/execution/RESULT.json. Verify the recorded artifact hashes. Reuse the
admitted interaction_decomposition_v1/actual_v1/RESULT.json, its exact literal
C/P/B coefficient dictionaries and paired reference targets, static evaluator
energies, and the original parameter domain. Pin all direct authorities before
execution. No outcome-based subset of origins or trials is allowed.

For each W8/W9 parent, C is the metal-shell interaction, P is the sum of literal
metal-water interactions, and B=C-P. Form Dy-minus-La separately at each shell
size. The four added errors are P8, B8, P9, B9, in kcal/mol. Only these paired
common-basis electronic-plus-original-D4 references have the required dispersion
cancellation. Absolute per-parent total-energy references remain unavailable.
Do not replace these targets with raw-cluster CC or uncorrected PBE0 differences.
C and the W9-minus-W8 differences are derived reports, not extra observations.

Retain the original 531-entry residual, four block scores and all 25 E/F gates.
Add descriptive J_P=(error_P8^2+error_P9^2)/2 and J_B analogously, using 1 kcal/mol
as an engineering scale only. Do not combine the six blocks into a new selection
objective. These are correlated development diagnostics, not a likelihood.

## Saved local response analysis

At every original origin, reconstruct derivatives from the exact saved plus/minus
vectors and their actual normalized parameter spans. Verify the old 531-row
Jacobian against the recorded one. Form the additional P/B derivatives from saved
energies and literal cycles. Retain every origin's gate status; do not describe
an origin with failed gates as qualified. The final selected vector differs from
the last Jacobian origin and must not be substituted for it.

For each residual block r with Jacobian J, compute g=(r^T J)/(r^T r), the
first-order fractional change of its residual norm. Solve three bounded linear
programs per origin: maximize t>=0 subject to g_P·dz<=-t and g_B·dz<=-t. Protect
(A) no other blocks, (B) the original E/F block, and (C) all four original blocks
by imposing g_protected·dz<=0. Use the original normalized parameter box and
the original +/-0.025 trust bounds. Use SciPy HiGHS with one thread, primal/dual
feasibility tolerance 1e-8, and retain primal feasibility, dual objective and gap
(required <=1e-7 for a numerical interpretation). No silent solver substitution.
If a block norm is <=1e-12, record the numerical degeneracy and do not assert
strict common descent from its normalized gradient.

Report the actual direction, each directional derivative, and the full linearized
residual norm at half and full proposed step. These directions are mathematical
probes only: they are not evaluated candidates. Nonincrease of an aggregate
E/F block does not preserve its 25 individual gates. A tangent solution does not
prove finite-step feasibility. A zero optimum is limited to this tangent problem
at this origin/domain; it cannot prove global failure of the functional form.

## Existing finite trials and numerical checks

For all 96 original trials, report actual versus linearized four component
errors and original residuals, including block-wise linearization discrepancies.
Record whether both J_P and J_B decrease, whether protected original blocks
increase, and whether all 25 stored E/F gates pass. Preserve failed-origin and
failed-trial distinctions. Numerical comparisons may use 1e-10 relative allowance
for equality; this is arithmetic handling, not a physical acceptance criterion.
Report all vectors descriptively, including individual errors and gate status.
Do not choose or promote an alternative historical iterate based on these data.

Verify the selected-vector decomposition reproduces the admitted result, literal
cycle closure and N-fold ion subtraction, vector/parameter joins and exact
derivative order/spans. Meaningful solver controls must include compatible and
opposed gradients and bound contacts. Independent review checks joins, arithmetic
and interpretation after execution. Preserve failures in distinct attempts.

## Resources and decisions

One numerical thread on the explicitly authorized login allowance; zero native
programs, zero new reference calls, zero scheduler allocation. Record wall and
CPU time, numerical library versions, source/input hashes and actual denominators.
Other user Slurm jobs are preserved. No infrastructure change is part of this run.

If existing trials or local responses show a credible joint improvement, the next
step is a separately declared finite nonlinear check of that mechanism. If the
responses expose a local conflict, use its specific observables to design an
additional interaction term or a different parameterization. Neither outcome
alone qualifies the model, establishes missing physics uniquely, or releases
fresh validation, bulk/protein sampling or an exchange prediction.
