# Fixed shared-pool comparison and reference operations

`scripts/nikasha_pool_compare.py` performs no molecular calculations. It checks
the algebra of actual collected matrices and preserves all source/status records.
The scientific selection rule and finite populations remain in [PLAN.md](PLAN.md).

## Operations

- `inspect --collections PATH [PATH ...] --output PATH [--reference PATH]`:
  report each real source's static, own-proposal and shared-pool contrasts, changes,
  selected candidate IDs and old-band transfer calls. Mathematical and operational
  pool variants remain separate. Without a compatible variant-specific reference,
  new decisions remain unavailable while computed raw contrasts remain visible.
- `calibrate --collections PATH PATH --agreement PATH --output PATH`:
  require the actual `pilot4` and `remaining26` result union to match the original30
  source IDs. Membership and labels come only from the released25 canonical context
  calibration record. Each variant independently uses the identical class-extrema
  rule and minimum gap strictly greater than0.02 model kcal/mol. Missing canonical
  results yield unavailable bands; no subset, crystals, PLM sources or folds enter
  calibration. Save this immutable artifact before pooled fold results are inspected.
- `compare --collections PATH --reference PATH --prior-comparison PATH --output PATH`:
  require the actual full225 source collection and its exact earlier proposal
  comparison. Retain all225 sources, the100 La- and125 Ca-conditioned arms, strict
  La4/Ca5/equal-arm medians and every100 previously declared La triple. Keep both pool
  variants under old and new bands alongside static context/native/core/composite,
  preserved DFT and old own-proposal results. Pairwise common-coverage counts and
  transitions are explicit. No reference fitting occurs on the fold data.

Each PATH is a required command argument, not a shell-dependent default. Outputs
must be new paths; the tools never overwrite a saved result or calibration. Actual
collection paths are supplied by the existing runner after completion. The prior
full fold comparison is already available at
`workspaces/nikasha_recovery_20260922/proposal_comparison.json`.

The new reference ID is `Nikasha_common_geometry_canonical25_reference_v1`, with
separate `_mathematical` and `_operational` variant IDs. It pins both source
collections, the original reference, software/model/method/settings, all25 values,
class extrema, any missing members and the freezing timestamp. Calibration accuracy
is not independent validation. Protocol IDs and production defaults remain unchanged.

## Actual checks and limits

The initial comparison suite uses the real eight-cell MACE-stage collection from
job1209840, with missing cross-solvent values left missing. Six tests pass and one
scientific-integration check is explicitly skipped until real cross-solvent results
exist. Existing canonical energies test only the unchanged extrema algebra; they
are never presented as new pool energies. A deliberately corrupted copy of the real
partial collection is rejected. Missing reference tests preserve the actual raw
score without replacing its calibration with zero or an inherited threshold.

Run the scoped suite with the existing driver:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python -m unittest discover -s tests -p test_nikasha_pool_compare.py -v
```

All populations are consumed structural repeats; unknown PLM labels remain unknown.
Ranges and fixed medians are structural descriptors, not equilibrium populations,
thermal entropy or quantitative affinities.
