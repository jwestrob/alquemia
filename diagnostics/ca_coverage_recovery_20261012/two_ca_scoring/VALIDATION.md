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

## Actual preparation and metadata binding

Job2969 subsequently completed: all five Q9 sources are prepared, each with192
local atoms, equal paired maps and the exact frozen canonical state key. Its
single-Ca regression passed. These are preparation results, not classifier results.

The actual preparation record is
`workspaces/ca_coverage_recovery_20261012/spectator_preparation_q9_v2/execution/FINAL_PREPARATION.json`.
Metadata-only binding then produced
`workspaces/ca_coverage_recovery_20261012/two_ca_scoring_q9_v1/manifest.json`
(SHA256 `41d8de1795fbac22fc661cc33361c1a8ec9a00ed97db3577ca7df9d6cf6bebf3`),
two lane manifests for samples0/2/4 and1/3, and `SOLVENT_CONFIG.json`
(SHA256 `25700e06ad665f696a48020fe20af262802f7268ffc70e801ee114d4df4f8667`).
No model or solvent call occurred during binding. Code, plan and wrappers are
frozen by those manifests. Root owns submission and actual integration results.

The native worker retains its existing two archived force/energy qualification
fixtures per lane. The same original kernel, optimizer, selector, model and
common-pool rules remain pinned. Corrected solvent collection retains existing
component-v2 acceptance and the failed historical absolute-v1 parity gate.
No production promotion, reference fit, new triple policy, or GPU extrapolation
is part of this implementation.
