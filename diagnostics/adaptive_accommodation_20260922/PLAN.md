# Archived-force active-coordinate diagnostic

Scope recorded before force projections, 2026-09-22. Authorized bounded parallel
task from root under Jacob's resumed accommodation goal. No molecular calls,
optimization, score modification, fitted labels or new reference.

Use every task in the immutable original30 `proposals_v1/manifest.json` and
primary225 `fold_proposals_v1/manifest.json` under
`workspaces/accommodation_nonlinear_20260920/`. Retain the 17 unprepared folds and
all failed/missing endpoint receipts. Analyze archived origin and final-proposal
force arrays only. Final endpoint geometries may differ; their residual loads
are separate diagnostics, not a gradient of Ca-minus-La at a common geometry.
Joint selection uses only matched origin coordinates. These are consumed sources.

For physical Jacobian J and mapped context Jacobian Jc, project the actual force
as g = -Jc:F, including source/cap chain rules. Report original radian/angstrom
derivatives. Normalize each mode by L = ||J_physical,heavy||_F, excluding synthetic
caps and H from this geometric metric. g/L has kcal/mol/angstrom units for one
angstrom of total physical heavy-atom displacement norm; it is not a relaxation
energy or stiffness. Do not divide by atom count or compare raw torque to force.

At common origins report signed d = gCa/L - gLa/L and each individual load.
A diagnostic four-mode preview uses supported angular primitives only: first
largest absolute differential load; then largest individual Ca/La residual load;
alternate these criteria for up to four slots. Modified Gram-Schmidt uses the
same physical heavy-atom metric, subtracting the same linear combinations from
both endpoint covectors; reject residual norm <= 1e-5 after unit normalization
(the existing 1e-10 squared-metric numerical rank scale). Deterministic mode-ID
ties. This is a proposed bounded selection rule, not a deployed model; no label
or score enters it. Report discarded and selected primitive IDs and residual
load coverage. Metal translations are an isotropic three-coordinate block,
reported separately, not three orientation-dependent competitors for one slot.

Also project final forces onto every existing mode at each endpoint's actual q;
report active versus omitted loads and archived boundary flags. Keep native
vacuum OMOL gradients distinct from composite gradients. For the already
completed seven-context analytic response, replay only existing mode derivatives
and normalize with its exact mappings to compare vacuum/composite rankings.
Do not infer solvent forces for the full reference folds.

Earlier broad-coordinate preparation already failed to improve native margins
and ended on all ten trust-region boundaries. This diagnostic must establish
what the restricted terminal modes leave unresolved, not rediscover that broad
vacuum minimization lowers energies. PQQ, water and added scaffold motions remain
explicitly unsupported until source-connectivity mappings and physical bounds
exist. No invented rigidity, curvature, entropy, populations or affinity score.
