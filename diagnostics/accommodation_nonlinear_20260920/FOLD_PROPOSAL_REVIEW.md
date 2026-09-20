# Primary225 prelaunch integration review

2026-09-20. Root approved the minimal shared-engine refactor before editing and
independently validated the resulting manifest before molecular execution:

`workspaces/accommodation_nonlinear_20260920/fold_proposals_v1/manifest.json`

SHA256 `0aedaee1064af04516224ff2cc6729843aa5a3e724a62b11477f6360eebb8cdb`.

Root reviewed the fixed source/mapping membership, 225 declared cases, 208 valid
original preparations, 416 paired tasks, unchanged scientific settings/reference
and CPU/GPU wrappers. The sole unavailable q0 is
`a8r3s4-pqq-la_model__conditioned_Ca__seed-1_sample-3__La`: its previously failed
native solvent calculation remains unavailable. The other 415 origins have exact
native-state, input and execution checks. No fresh q0 is scheduled.

Root's explicit coordination decision: after the seven-test rerun and independent
artifact review pass, submit this GPU manifest and then all four actual solvent
shards after their dry-runs. No further permission round is required. This implements
Jacob's standing overnight research authorization under FOLD_TRANSFER_PLAN.md;
it does not promote a production method.

The first seven-test attempt retained a metadata-only assertion failure: original30
collection adds declared source-group metadata. Every original scientific field
replays exactly. The corrected regression requires those old fields exactly and
restricts added keys to root_case_id, biological_group, source_conditioning_metal
and canonical_coordinate_match. The first passing rerun is recorded in
FOLD_PROPOSAL_TESTS_v2.txt; the additional metadata-key restriction is checked in v3.
The failed v1 test receipt is preserved. No scientific inputs or implementation
snapshot changed in response to this assertion.

The shared engine's default30/60 policy, optimizer, gradients, energy selection
and original result values are unchanged. Only explicit primary225 cardinality,
four-shard staging/receipt lookup and reporting metadata were added. Each CPU
shard has separate contained paths/lock; the master has zero executable tasks.
Missing proposals never borrow the archived baseline value. The separate comparison
retains those archived values in their own unchanged fields.
