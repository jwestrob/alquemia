# Independent review: finite intact-PQQ pose test

**No scientific launch blocker in the declared finite plan.** This is a fair
exploratory test of one omitted cofactor coordinate, with a useful stopping rule.
It is not a comprehensive test of PQQ relaxation or a diagnosis of all nine
inconclusive sources. Review is read-only; no calculations or submissions.

## Why this matrix is informative

The saved force result4598 separates a largely shared carbonyl-shortening load
from a more source-dependent rigid e3 rotational response. The saved Ca−La e3
projections have different signs between samples0/3, both at origin and along
adaptive_La. These are MACE generalized forces, not composite gradients or
independent DFT forces. They motivate measuring a finite response; they do not
predict its curvature or a desired class direction.

The declared ±1° rotation moves intact27-atom PQQ about its24-heavy-atom centroid.
It preserves chemical state and internal geometry. Both signs, both metals and
both sources are retained. Applying it to **all three** existing donor geometries
makes24 new endpoints and12 reused baselines. Comparing complete common3→9 pools
avoids calling a favorable single path useful when another original geometry
already wins. Adding adaptive_Ca, whose force was not in the eight-cell projection
study, completes the existing pool rather than choosing a favorable subset.

Preparation must retain the specified exact paired coordinates, unmoved exterior
and caps, internal/H-parent distance invariance, inherited contact gates and actual
maximum displacements. Missing required cells must keep the full comparison
unavailable. The plan states these conditions; this review does not certify their
future numerical implementation or outcome.

## Interpretation and stopping rule

- **Large change:** demonstrates that this finite pose coordinate changes the
  within-model paired vacuum minimum on the fixed donor pool. The matched solvent
  stage is then justified. It does not establish that the original pose is wrong,
  that the alternate pose is an equilibrium minimum, or that the PQQ force is
  physically correct. Prior ExaA DFT validated donor work, not this new PQQ motion.
- **Small change:** supports closing this particular ±1° e3 vacuum-pool probe.
  It does not exclude other axes, different angular scales, coupled donor motion,
  internal PQQ strain or a solvent-dominated pose effect. Because solvent is
  intentionally screened using vacuum response, an unexecuted solvent stage must
  be reported as untested, not inferred small.
- **Large common change:** is not source-discrimination repair. Preserve each
  metal's work and selected geometry, each source's paired shift and the
  target-minus-control change. More favorable individual energies need not
  improve the paired contrast; extending a pool can lower both rows similarly.
- **Continuation threshold:**1kcal/mol inherits the prior declared interpretive
  resolution. It is a pragmatic resource gate, not uncertainty, a numerical
  admission tolerance or a physically derived cutoff. Either source and either
  sign can trigger solvent. No label, expected class, band crossing or fitted
  threshold selects continuation. This is appropriate for the bounded test.

No further fractions, optimizer development or calibration change follows
implicitly. A substantial matched composite response would still need independent
physical support and the designated reference/transfer population before any
classifier claim. This branch and the independently referenced MMOL motion answer
different questions.

## Evidence reviewed

- [Frozen pose plan](../PQQ_POSE_PLAN.md).
- [Actual saved-force interpretation](../pqq_force_inventory/RESULT_REPORT.md) and
  [copied numerical summary](../pqq_force_inventory/RESULT_SUMMARY.json).
- Underlying saved-force result SHA256:
  `0546dd64057effee7c13ad01c7d7541a73c6b9bf3776b38dad4e8e08e2df7eb6`.
- Saved-force execution SHA256:
  `9e1482bfd16b3def89f533d4bc60cbd1c82b3d6df2dfba9d9c93d94903513e82`.
