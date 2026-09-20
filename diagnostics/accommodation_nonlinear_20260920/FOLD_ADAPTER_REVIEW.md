# Independent primary225 adapter review

2026-09-20, water_basins agent. **No blocking issue found.** This is an
implementation/input audit, not a new scientific result or launch receipt.
No molecular calculations or source changes were made by this review.

Reviewed `accommodation_fold_proposals.py`, the shared proposal-engine diff,
the frozen transfer plan/reference, and the actual prepared manifest:
`workspaces/accommodation_nonlinear_20260920/fold_proposals_v1/manifest.json`,
SHA-256 `0aedaee1064af04516224ff2cc6729843aa5a3e724a62b11477f6360eebb8cdb`.

## Actual artifact checks

- All 225 original primary source IDs remain: 125 Ca-conditioned and 100
  La-conditioned; 208 supported preparations plus all 17 unavailable records.
  No canonical calibration replay entered the transfer population.
- The 416 prepared endpoint tasks retain original XYZ pins, charges,
  multiplicity and common Ca/La physical mappings/active coordinates.
- Exactly 415 q0 endpoints are available. The original failed
  `a8r3s4-pqq-la_model__conditioned_Ca__seed-1_sample-3__La` remains null and
  prevents that endpoint's proposal. It is not replaced by a baseline success.
- Independently checked 830 retained q0 low-level input recipes and parsed
  coordinates. Their receipts point only to executed
  `solvent_shards_v2/shard_0..3`, not the unexecuted parent manifest.
  Every one of the 207 complete archived pair contrasts replays exactly.
  The expensive original output parsing was not repeated by this review.
- An AST comparison against commit `bdec940` finds `MACEProposal`, `optimize`,
  `select_endpoint` and `score` unchanged. Bounds, objective, starting point,
  optimizer and selection decrease remain the same as the 30-source experiment.
- Old bands and `PQQ_MACE_donor_proposal_composite_canonical_reference_v1` are
  separate pinned artifacts. The latter uses only the designated canonical25;
  no transfer-fold score enters calibration. Strict La4/Ca5 and all100 triple
  summaries preserve missing members and do not choose the best subset.

## Runner and reporting

Four solver shards place each task inside its owning manifest directory, giving
each executor separate locks/events. Task-to-manifest pointers resolve actual
receipts through that shard. The parent manifest has no executable tasks.
Actual displaced-input manifests are generated only after the GPU proposals;
their mandatory existing-runner dry-runs must still pass before submission.

One minor reporting issue was fixed during review: archived-origin replay fields
now read the fold comparison for this scope instead of searching the canonical
reference table. This does not alter selection or numerical scores. The separate
comparison retains original DFT/native/composite scores even when the proposal is
unavailable. Failed or missing proposal energies never become successful origin
fallbacks; origin selection requires both actually evaluated geometries.

Seven real-fixture tests pass, zero skips, in 35.850 s:
`FOLD_PROPOSAL_TESTS_v3.txt`. They include real-input four-shard staging/dry-run,
archived shard receipt lookup, exact pair algebra, complete source membership,
prelaunch missing-result accounting and unchanged 30-source selection fields.
The initial regression failure concerned newly added biological-group metadata
only; the corrected test requires every original field to match exactly and
permits only the four explicitly declared new metadata fields.

This audit supports executing the already authorized finite manifest. It does
not establish better discrimination, a native minimum, or production promotion.
