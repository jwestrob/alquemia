# Static masked-MACE context candidate rejected before new inference

The exact cached GGR partition contrast is20.160308kcal-scale, failing the
predeclared2kcal-scale tolerance. **Do not launch the six conditional whole-protein
calls for this candidate.** The production baseline and prior scores are unchanged.

The model was DFT/CPCM(core) + masked-MACE(full) - masked-MACE(core), with fixed
unit coefficient and no new solvent term. For connected-minus-extended GGR,
the identical full protein cancels algebraically:

| Component | Partition shift |
|---|---:|
| Native DFT core R | -7.343500873 |
| Masked learned core R | -27.503809106 |
| Hybrid R: DFT shift minus learned shift | +20.160308232 |

DFT quantities use kcal/mol; learned quantities use the explicitly converted
model-kcal scale. Unrounded endpoint energies and both terms remain in the
[result](MASKED_SUBTRACTIVE_CONTEXT_RESULT.json). The raw-zero feature does not
supply a stable core-replacement descriptor under this partition test. No full
protein term can repair this particular discrepancy because it is common to
both representations. This rejects this fixed candidate, not every hybrid model.

All four declared core representations and their matched original-H full systems
pass source-atom, element, cap and paired-coordinate checks; maximum residual
is4.99e-11A. Both GGR representations point to the same whole-system coordinates.
Known original-source H defects remain a limitation. No H-normalized full-chain
score was silently combined with these cores. There was no new DFT, MACE,
solvent or training call; local receipt replay/analysis time is recorded separately.
Two actual-fixture regression tests pass in4.384s, including direct endpoint
algebra and rejection of a corrupted copy of a real score report.

The same learned descriptor also failed the preceding unanchored response and
curvature screens. Neither local relaxation nor this static subtraction is
promoted. Its earlier PQQ calibration and failed GGR robustness results remain
unaltered. The goal continues with a separately declared [shared learned neutral
feature](SHARED_NEUTRAL_FEATURE_PLAN.md), one fixed conditioning control with no
charge-category sweep or fabricated physical neutralization.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_omol_context.py --response-report workspaces/mace_omol_20260917/masked_response_report_v1/result.json --agreement diagnostics/mace_omol_20260917/MASKED_SUBTRACTIVE_CONTEXT_PLAN.md --output workspaces/mace_omol_20260917/masked_context_partition_review_v1
```

This command verifies and replays archived results only. Use a new output directory.
