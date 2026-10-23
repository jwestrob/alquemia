# Saved local response diagnosis

## Question

The completed four-parameter pilot17567 reached its coordinate step stopping rule, with La and paired maximum-component errors still above5 kcal/mol/Å. Does its already-sampled local parameter response predict a useful coupled step, or does the remaining weighted residual mostly lie outside that local response? This is a saved-data diagnosis, not another fit or a parameter-family verdict.

## Fixed evidence and arithmetic

Bind the terminal pilot result SHA256 `524d093fb920304e254369097a96843408ab21597aa7898f5d1a3947c64fe6b6`, original SPEC/objective/search and four-parameter family. Preserve its selected candidate exactly. Use the eight existing final sensitivity vectors only to construct a central/clipped finite-difference Jacobian of the **same weighted residual vector** with respect to the four normalized parameter coordinates. Denominators are the actual saved normalized displacements; no new evaluations or tolerance changes.

Use SVD least squares with relative singular-value cutoff1e−6, matching the pilot's identification convention. Report J singular values/rank, the unrestricted Gauss–Newton solution to `min ||r+J delta||²`, and the projected residual split into the part within the column space and its orthogonal remainder. Report loss and force/energy block decomposition; this linearized remainder is not a lower bound on the nonlinear model family.

Also report two bounded predictions: uniformly scale the unrestricted direction so no normalized component exceeds0.025 and it stays inside the original parameter box; and solve the local linear least-squares problem directly within `max|delta|≤0.025` intersected with the original box. The latter uses exhaustive four-coordinate active-set enumeration (81 free/lower/upper patterns), solving each free least-squares subproblem, admitting feasible solutions with machine-only1e−12 slack and choosing the minimum predicted loss in fixed enumeration order. These are hypothetical predictions, not evaluated or selected candidates. Report normalized and physical parameter displacements separately.

Construct the analogous finite differences of saved raw body-response residuals so every predicted step reports the original separate La/Dy/paired RMS and maximum-component checks plus balanced energy error. Keep the same block weights, geometry/category averages, force and torque scales. Verify their weighted representation against the saved residual vector at machine precision (absolute1e−10 plus relative1e−10); this checks bookkeeping, not molecular accuracy.

For every already-evaluated vector, retain its distance from the selected point, predicted and actual loss and separate controls. Mark the prespecified local neighborhood `max|delta|≤0.025` and the eight sensitivity points. Explicitly show prediction error there; far-away points do not validate a local linear approximation. Do not replace the selected candidate even if a previously unused sensitivity point has a smaller actual objective.

## Scope and decision limits

Run saved arithmetic on one local CPU thread using NumPy linear algebra, with all numerical thread controls1. Record source/input hashes, package version, wall/CPU costs and zero native calls. No native engine imports, new geometries, epsilon changes, altered bounds/weights, validation reads, refit, retry or production release.

A predicted coupled improvement motivates, at most, a separately declared finite native test after root review. Poor local predictions, weak response or an orthogonal remainder can identify a local limitation; none proves global nonrepresentability or licenses a larger family automatically. The uncorrected PBE0 target and CC/CP/state limitations remain unchanged.
