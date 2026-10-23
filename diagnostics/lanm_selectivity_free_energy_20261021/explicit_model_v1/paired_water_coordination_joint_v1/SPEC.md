# Joint coordination development implementation specification

The immutable scientific declaration is ROOT_DECLARATION.md. This implementation
extends the prior eight-cell water fit to exactly ten consumed cells, adding the
original W9 La/Dy configurations and keeping both old and new individual checks.

## Data and computation

The evaluator prepares a ten-cell input/reference authority and replays the frozen
17605 selected eight cells plus the exact17670 W9 outputs. Preparation performs
zero native calls. Parameters, bounds and numerical engines retain their original
six-parameter definitions. Effective input reuse depends on the represented metal,
three water-pair values, geometry, exact key/parameter bytes and native source.
No orientation or protein outcome enters the objective.

The driver preserves the old six-dimensional derivative/proposal algorithm and
uses the joint objective. The old eight-cell score is independently reproduced
at the initial point and reported separately at the final point. The complete
joint inventory has413 residual components,10 energy controls and15 force controls.
Energy/force acceptance is per control; lowering the average is insufficient.

## Finite execution

One task, eight CPUs, eight single-thread workers, standard_expand and scheduler
default memory/time. Maximum24 outer steps,384 new full parameter vectors and
5120 new analyze/testgrad program invocations. Initial ten-cell replay is free of
native calls. A derivative batch changes only one represented-metal triple,
allowing unchanged metal cells to be reused. Four trial proposals may change both.

Derivative probes cannot be selected. Evaluate the fixed trial batch completely;
strictly lower objective wins, and current point wins ties. Stop on all individual
checks passing, no improvement after the complete batch, native failure, or the
finite cap. No automatic second start or extension. Report exact native counts,
selected bounds, stop reason and the separate origin of the last Jacobian.
NumPy rcond=1e-6 is a relative singular-value cutoff.

## Admission and interpretation

Actual zero-call preparation, saved objective controls and independent review
precede root submission. Preserve release-pinned pre-submit records permanently;
write post-run observations to distinct new filenames. Scheduler completion and
scientific acceptance are separate. Failure remains data; it does not prove global
nonrepresentability. Success only establishes this consumed-set fit, requiring
fresh transfer testing and matched donor reevaluation before further promotion.
