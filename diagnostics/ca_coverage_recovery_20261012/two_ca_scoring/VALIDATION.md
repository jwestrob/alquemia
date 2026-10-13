# Implementation status

The native and corrected-solvent adapters are implemented. Python AST parsing
passes for both scripts; `bash -n` passes independently for all three Slurm
wrappers. These are syntax checks only. No scientific import, molecular parser,
model call, geometry evaluation, or integration test ran on the login node.

Root reviewed the code. The preparation interface is matched to Curie's actual
schema: canonical anchors come from `PREPARATION.cases`, selected by both
`anchor_case_id` and `selection_id`; admitted rows require a passed spectator
receipt and exact frozen state key. The real single-Ca regression receipt must
pass before an eligible batch can be bound.

Q9 preparation job2969 is pending/running at this documentation checkpoint;
actual admission and scoring results are not yet available. No runnable source
manifest is inferred from successful folding. The first admitted batch uses two
independent seven-CPU lanes. The other ten samples use a separate later batch
with this same implementation. Scientific integration remains **UNRUN** until
actual Slurm receipts demonstrate otherwise.

The native worker retains its existing two archived force/energy qualification
fixtures per lane. The same original kernel, optimizer, selector, model and
common-pool rules remain pinned. Corrected solvent collection retains existing
component-v2 acceptance and the failed historical absolute-v1 parity gate.
No production promotion, reference fit, new triple policy, or GPU extrapolation
is part of this implementation.
