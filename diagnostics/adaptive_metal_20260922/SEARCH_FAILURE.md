# MMOL1770 La: a rejected search trial, not unscorable source chemistry

Read-only diagnosis of unchanged angular-continuation job1209886, task
`mmol_1770-pqq-la_model__La`. Its original failure remains untouched.

The source q0 maps correctly, returns a successful native-MACE energy, and
replays the archived origin exactly. The next SLSQP callback repeats q0 from
cache. The third request—the first nonzero proposed step—is
`[-0.8, 0.8, 0.30437587274943945, -0.8]` rad on
`A/318/chi2, A/362/chi1, A/200/chi3, A/200/chi1`.

That step moves a source heavy atom **2.525832289949 Å** and brings context
indices **21/85**, actual **Glu200 OE1 / Trp322 HE1** atoms, from
2.597638350502 Å to **0.547794072825 Å**. The unchanged H-containing-pair guard
requires avoiding a newly introduced separation below **0.55 Å**. Thus it
rejects this trial before a MACE request. Bond preservation, fixed atoms and
analytic-map checks themselves pass. No iteration is accepted and no candidate
energy exists for that rejected geometry.

The spherical heavy-atom inequalities are `0.8²−|Δx_i(q)|²≥0`. Their Jacobians
are exactly zero at q0, so the first linearized SLSQP subproblem does not feel
the finite sphere. Although the approved algorithm permits infeasible trial
coordinates, its safety guard raises an exception instead of returning a
numerical function value. That exception ends the search before SLSQP can
contract its step. It is a search/interface failure; the valid source energy
and a feasible neighborhood remain available.

## Geometry-only contraction of the actual failed direction

| Fraction of failed step | Maximum source-heavy displacement / Å | Overlap check | Within 0.8 Å |
|---|---:|---|---|
| 1 | 2.525832290 | fails | no |
| 1/2 | 1.314484745 | passes | no |
| 1/4 | 0.667086339 | passes | yes |

These are deterministic coordinate checks only. No displaced energy was
evaluated, and a feasible point does not imply that it lowers energy or provides
a converged solution. Exact source/failure pins and checks are preserved in
`workspaces/adaptive_metal_20260922/prepared_v1/MMOL_FAILED_STEP_GEOMETRY.json`.

## Proposed next numerical change — not implemented

Use a driver whose **accepted-state line search owns the trial coordinates**:
retain the last feasible accepted q, generate a descent/SQP direction using its
actual energy/gradient, and halve that direction until the exact physical map,
component bounds and overlap guards pass. Only then request a real native energy
and use actual energy decrease in line acceptance. Retain rejected geometric
trials and declare a search failure if no usable step exists. The same source,
physical domain and objective remain; no proxy energy is needed. A version using
this approach must define its constrained direction and boundary/stalling rule
before molecular execution and test them on the archived failed step.

**Do not implement contraction inside the current SLSQP objective callback and
return the contracted point's energy/gradient as if evaluated at SLSQP's requested
q.** Those quantities would describe different coordinates and break the
optimizer's mathematical contract. Similarly, do not make up a repulsive energy
for the rejected point, loosen the collision threshold, clip a final candidate
into success, or silently restart the original continuation.

The current joint-metal adapter intentionally preserves the old numerical
behavior so root can review this change separately. The first-stage numerical
question is handling rejected trial geometry coherently; it does not require
changing chemistry, fitting labels or discarding the unavailable denominator.
