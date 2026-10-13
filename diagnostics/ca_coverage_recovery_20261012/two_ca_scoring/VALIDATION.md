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

The remaining ten actual C5AXV8/Q88JH5 preparations also passed and were bound
without changing the adapters: `two_ca_scoring_c5q88_v1/manifest.json` has SHA256
`38c4630f8645ab31bbd88870b3d101f89b22cb40fe8215bbd10f33e40ca36ea2`;
its three seven-CPU lanes contain4/3/3 sources. Its solvent configuration has
SHA256 `0824aad4ca5f7ecf1b4ff52e091148fe9ece305b4bed36577efe6300b9bdb63f`.
Root submitted Q9 lanes2996/2997, collector2998 and solvent2999; remaining
lanes3018/3019/3020, collector3021 and solvent3022. Read actual receipts for their
current states rather than treating submission as scientific acceptance.

## Cumulative collection

Separate new files `MERGE.py` / `merge.sbatch` leave the running code and pinned
plan untouched. AST/shell checks pass. The metadata-only final manifest is
`workspaces/ca_coverage_recovery_20261012/two_ca_scoring_final_v1/manifest.json`
(SHA256 `17fcbd7bf20429e468cba0aaddf72d354c6e2acb1e63feab9e8b88a4fd0fcf10`).
Root can run it on one Slurm CPU after both solvent jobs terminate. It requires
both real5+10 result overlays, validates their exact original125 and prior-three
records, and writes the cumulative JSON/TSV/report. A missing batch explicitly
fails with an execution receipt; it cannot be reported as a completed15-source
comparison. This collection has not yet executed at this checkpoint.

The native worker retains its existing two archived force/energy qualification
fixtures per lane. The same original kernel, optimizer, selector, model and
common-pool rules remain pinned. Corrected solvent collection retains existing
component-v2 acceptance and the failed historical absolute-v1 parity gate.
No production promotion, reference fit, new triple policy, or GPU extrapolation
is part of this implementation.

## Actual completion — 12 October 2026

The pending/unrun statements above are historical. All five native lanes and
both solvent batches completed; each native lane passed its two archived
energy/force qualification fixtures, and all required component-v2 gates passed.
The cumulative merger3049 executed successfully on Slurm. Independent review
confirmed exact primary125 preservation, the three prior recoveries and disjoint
five-plus-ten additions; JSON/TSV decisions and all scheduler rows agree. See
../DELIVERY.md for results, failures, costs and physical limitations. No extra
molecular computation was run to write this completion note.
