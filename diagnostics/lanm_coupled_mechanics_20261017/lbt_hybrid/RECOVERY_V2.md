# Derived cap-Jacobian replay recovery

## Actual cause

Frozen evaluator a482559c passed CPU preflight 12972 but failed scout 12975
before molecular calls on the H200 host. Exact-host zero-call diagnostic
12984 compared all four origin maps twice. Three replayed exactly. Only
7CCN/baseline differed: 32 duplicated scalar leaves in cap 0/core row 170's
retained and omitted Jacobians, at most 4 binary64 ULPs and
1.1102230246251565e-16 absolute. All other fields, including physical/core/field
coordinates, charges, identities and types, were exact. Repeated maps matched.
The original failure, evaluator and preparation remain unchanged.

## Narrow versioned change

EVALUATOR_V2 uses MAPPING_REPLAY_V2 for origin and restored-point mapping
replay. The helper requires identical tree structure, types, cap/core duplicate
rows, identities, charges and coordinates. Only float coefficients of the two
3×3 retained/omitted analytic cap Jacobians may differ, at both duplicated
locations. Each difference must be at most 8 binary64 ULPs **and**
8 machine epsilons absolute (1.7763568394002505e-15). Nonfinite values, changed
zero values/signs, changed dimensions and any other leaf difference fail.

This is a small bounded floating-point recomputation allowance for dimensionless
analytic coefficients, with a factor-of-two ULP margin over the observed
four-step difference. It is not a molecular-accuracy threshold, a physical
force tolerance or a guarantee on every host. Any larger difference still
stops execution. No comparison of energy outcomes informed this choice.

The archived map remains the evidence for an existing point. Restored component
gradients are recomputed with that archived map and retain exact raw replay.
New component calculations continue to use the original unmodified mapper at
their actual coordinates. The exact task-XYZ check is unchanged. All physical,
state, energy, SCC, cache-coordinate and finite-difference gates remain unchanged.
The recovery helper, diagnostic evidence and recovery specification enter the
new cache identity. No prior V1 molecular results exist to migrate.

## Required zero-call preflight

Before the nine-composite scout, root runs V2 preflight on the actual failed
host with the same Python and launcher environment. It must replay all four
actual diagnosed maps and independently regenerate/replay all four current
origin maps, retaining every allowed coefficient residual. Existing eight
native GFN/solvent replays and frozen finite-difference directions still run.

Twelve negative controls must reject: one-ULP changes to physical, field,
real-QM and cap coordinates; one-ULP external charge change; physical/core
identity changes; Jacobian float-type change; a nine-ULP Jacobian change;
1e-8 Jacobian change; nonfinite Jacobian; inconsistent duplicated cap row.
Cap-coordinate and Jacobian controls modify both duplicates together except
the dedicated inconsistency control, so rejection tests the intended gate.
These are saved-data tests only, with zero model or molecular calls.

The unchanged finite nine-cell scout follows only after actual preflight
acceptance. Root owns all submissions; no production change is made.
