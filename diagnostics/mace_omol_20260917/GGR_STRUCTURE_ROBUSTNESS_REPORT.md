# Masked MACE fails the expanded GGR robustness test

All four new forwards completed and passed numerical accounting. The frozen alpha-minus-GGR ordering passes only **2/6 structural comparisons**; the all-case robustness gate is **false**. These are two biological groups, not six independent affinity tests. The earlier five-case development pass and the separate PQQ calibration remain valid historical results.

| Structure | Descriptor, model kcal |
|---|---:|
| GGR_1GLG | 23.987421458 |
| GGR_2FW0 | 44.884127357 |
| GGR_2FVY | 45.975500410 |
| ALPHA_1F6S | 36.729109306 |
| ALPHA_6IP9 | 28.099819820 |

GGR spans 21.988078952 model kcal. Both new structures exceed both alpha scores, reversing the desired ordering. No PQQ band, affinity zero, changed label or geometry rescue was applied.

## Diagnosis and limits

All three GGR inputs retain the same six donor residues and seven coordinating oxygen atoms. Several metal–oxygen distances differ by about 0.1–0.2 Angstrom. Saved readouts place most of the score change in the metal and nearby donor residues; the 6–12 Angstrom shell also changes appreciably. Readout closure passes for all five structures. These learned atomic contributions are bookkeeping, not observable energies or proof that one particular donor causes the failure.

The sources differ in sugar-associated conformation, terminal coverage and preparation. The new structures omit four terminal residues, including two lysines: physical protein charge is −8 versus −6 for 1GLG; those terminal omissions are over 28 Angstrom from the metal. Total-charge masking and vanishing remote readout differences do not establish a generally correct electrostatic response. No unique physical cause has been identified.

## Implementation and actual cost

Added a source-backed preparation bridge with explicit raw-structure inventory, exclusions, metal identity and peptide connectivity checks. It reproduces original 1GLG physical coordinates exactly and rejects raw 1KB0 TRO512 before filtering. New preparation policy: `omol_source_backed_chain_A_ff19sb_H_explicit_exclusions_v1`. The descriptor protocol and baseline remain unchanged.

Jobs 1200845/1200846 each used one A5000 for 57 seconds: **114 GPU allocation-seconds**, 1824 allocated core-seconds, 118.469 reported actual CPU-seconds. Four model calls took 52.124891 seconds total; peak allocated GPU memory was 6,108,277,248 bytes. Local preparation, tests and reports have separate receipts. No new DFT, training, solver, forces or relaxation calculations.

All inputs, unrounded energies, native readouts, attempts and scheduler receipts are pinned in the compact result and full workspace report. The comparison can be rerun without inference using the command below.

```bash
MACE_DRIVER=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
MACE_WORK="$PWD/workspaces/mace_omol_20260917"
"$MACE_DRIVER" scripts/mace_omol_ggr_robustness.py \
  --development-report "$MACE_WORK/charge_ablation_report_v1/result.json" \
  --ggr-2fw0-report "$MACE_WORK/ggr_masked_2fw0_report_v1/result.json" \
  --ggr-2fvy-report "$MACE_WORK/ggr_masked_2fvy_report_v1/result.json" \
  --agreement diagnostics/mace_omol_20260917/GGR_STRUCTURE_ROBUSTNESS_PLAN.md \
  --diagnosis-plan diagnostics/mace_omol_20260917/GGR_SAVED_READOUT_PLAN.md \
  --sacct "$MACE_WORK/ggr_structure_sacct_v1.tsv" \
  --output "$MACE_WORK/ggr_structure_review_v2"
```

## Decision

Retain the production baseline. The masked MACE descriptor is numerically credible and inexpensive on these inputs; broad affinity usefulness remains unestablished, and its apparent non-PQQ improvement is not structurally robust. Keep it available as a research PQQ classifier. Continue development without declaring this failed expanded test a success.
