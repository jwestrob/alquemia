# Numerical core for the one declared component campaign

ROOT_DECLARATION.json is authoritative. CORE.py has no file, native, scheduler or
reference API. The driver verifies immutable inputs, native admissions, budgets,
release, and cache identities. Only `propose()` calls SLSQP, exactly once after
validating its current point and24 probes. Preparation controls replace that
solver with a mock; no real campaign proposal is computed before release.

## Interface

- `make_fixed(origin, manifest, decomposition, domain, declaration)` validates
  the exact origin24,46 variable/48 static terms, literal paired C/P/B cycles,
  original-domain limits, fixed component denominators and four old ceilings.
- `assess(vector, fixed)` reports all25 physical gates, four ceiling flags, four
  signed component errors, J_P/J_B, q and the declared full development target.
  It does not certify native execution. The separate substantive-recovery flag
  requires q<=0.25, both component RMS values below the admitted original initial
  model, and all25/four ceilings. The quarter-q flag alone is not that milestone.
- `probe_points(parameters, domain)` returns the original ordered minus/plus24
  points, +/-0.005 normalized and clamped only in probe construction.
- `linearize(current, probes24, fixed)` checks exact probe identities and actual
  nonzero spans. It reconstructs the original413/531 rows, physical energy/body
  gates, and paired-component secants from the complete saved vectors.
- `propose(current, probes24, fixed)` starts at d=0,q=q_current, calls SLSQP once
  with analytic derivatives, maxiter200/ftol1e-12, and reports the original solver
  status plus four ordered fraction certificates. It never executes a vector.
- `select(current, trial_slots, fixed)` accepts the same four slots after the
  driver adds complete admitted vectors to eligible entries. Every eligible
  slot must be complete and explicitly native_admitted; otherwise selection
  fails. Certified parameters/key must match exactly. Probes are not candidates.

A complete vector dictionary has `status`, `parameters`, `model`, `response`,
`score` and all46 `cells`, following the frozen scientific evaluator. The core
keeps authoritative stored score/chi values and checks their arithmetic, avoiding
replay through the historical exact derived-LAPACK cache branch.

## Fixed constraints and independent certificates

Minimize q subject to each linearized P/B squared residual norm divided by its
fixed origin24 J_P/J_B being at most q. Each old residual block retains its fixed
origin24 ceiling, including the exact half-weighted531 slices. Retain ten |energy|
limits1 and fifteen force RMS/max limits2/5 with unchanged category weights.
These are851 affine/convex quadratic inequalities plus the original box/trust
bounds and0<=q<=q_current. Ceilings never become per-step monotonic constraints.

Each fraction f in1,.5,.25,.125 uses d_f=f*d_returned and
q_f=(1-f)*q_current+f*q_returned. A separate scalar implementation uses math.fsum
and the raw physical residuals; it does not call the vectorized constraint routine.
It verifies exact finite parameter/normalized trust bounds, all851 inequalities,
all25 predicted gates, old ceilings and strict q improvement. A successful solver
still needs this certificate. Failed solver statuses remain failed and require
strictly positive residual-constraint margins for a fraction to qualify; parameter
bound equality is allowed. No physical-threshold relaxation, alternate fraction,
solver retry or clipping of returned candidates is performed.

The actual selection epsilon is fixed at1e-10*max(1,abs(q_current)). Scan eligible
complete candidates in declared fraction order. Replace the incumbent only for
a lower q beyond that epsilon and all25+four fixed-ceiling passage. Current wins
ties, then the earlier fraction. One component or old block may worsen relative
to the previous iterate while remaining below its fixed origin ceiling.

No qualified fraction/no actual improvement/native failure triggers the driver's
frozen stop rules. Neither finite failure nor SLSQP status establishes global
infeasibility. A local witness is not a nonlinear physical-accuracy result. No
fresh validation, new reference, fit extension or production release is implied.
