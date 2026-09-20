# A distinct reference for the geometry-selection descriptor

The completed MACE proposal experiment preserves class ordering for all original
PQQ references. Its scores are a new descriptor: the old calcium band leaves
Q9Z4J7 inconclusive. The original frozen-band result remains in the proposal
report. It is not changed retroactively.

Before running the new fold-transfer experiment, a separate developmental
reference was frozen using the same 25 designated canonical calibration cases
and the original class-extrema rule. The three crystals and two PLM predictions
are excluded. No noncanonical fold scores enter calibration.

- Reference: `PQQ_MACE_donor_proposal_composite_canonical_reference_v1`.
- Protocol: `native_MACE_proposal_primary_composite_selection_v1`.
- Ca-supported: R <= -405462.0237183686 model-kcal/mol.
- La-supported: R >= -405456.1793213965 model-kcal/mol.
- In between: inconclusive. Gap 5.8443969720974565 model-kcal/mol.
- Same existing minimum calibration gap: 0.02 model-kcal/mol.

The 25 calibration cases are separated by construction and are not a validation
result. All three consumed crystals retain their expected calls. The unknown PLM
predictions remain one calcium-supported and one inconclusive under either set of
bands. A favorable shift or transferred label alone cannot establish improved
specificity prediction. Their true classes remain unknown. The known-class fold
challenge is the next utility test.

Full canonical inputs, actual selected geometries, energy-selection receipts,
implementation hashes, original bands and both sets of initial decisions are in
`workspaces/accommodation_nonlinear_20260920/proposal_calibration_v1/REFERENCE.json`.
The existing scorer, references and production results remain unchanged.

Reproduce into a **new** output path:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
diagnostics/accommodation_nonlinear_20260920/calibrate_proposals.py \
  --result workspaces/accommodation_nonlinear_20260920/proposals_v1/final_1204171.json \
  --agreement diagnostics/accommodation_nonlinear_20260920/FOLD_TRANSFER_PLAN.md \
  --output workspaces/accommodation_nonlinear_20260920/proposal_calibration_replay_v1/REFERENCE.json
```

Four tests pass on the real completed artifacts: component/selection replay,
calibration unchanged when transfer records are removed, missing-member status,
and rejection of corrupted labels/selection/identities. These are parser/algebra
tests, not fresh molecular calculations. No new energy calls were made here.
