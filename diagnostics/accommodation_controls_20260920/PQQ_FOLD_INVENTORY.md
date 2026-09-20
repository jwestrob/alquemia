# Existing canonical PQQ fold samples — 2026-09-20

**All 250 real archived AF3 samples are frozen for structural-robustness testing. No new folding is needed.** Primary evaluation pools contain 100 unselected La-conditioned samples and 125 Ca-conditioned samples; the 25 canonical La geometries are retained for replay and excluded from those pools. Every protein has seed 1 samples 0–4 under each conditioning arm.

## Ready files and preparation gate

- [All 250 fast-preparer source configs](PQQ_ALL250_SOURCES.json): exact existing config, residue roles, group and labels; unique case IDs; `archived_*` comparison fields retained only for the 25 canonical replays.
- Metadata: `workspaces/accommodation_controls_20260920/pqq_all250_metadata.json` — original reference, conditioning, sample, canonical match, primary-pool membership, confidence pins and readiness.
- Full inventory: `workspaces/accommodation_controls_20260920/pqq_fold_samples_inventory.json`; flat table `pqq_fold_samples.tsv`.
- Parent's frozen evaluation plan: [FOLD_ROBUSTNESS_PLAN.md](../accommodation_goal_20260920/FOLD_ROBUSTNESS_PLAN.md).

**125 La-conditioned sources are compatible with the current normalizer. The 125 Ca-conditioned sources require an explicit residue-name adapter.** All 250 pass the existing complete-PQQ source detector. No protonation, carving or scoring was performed, so downstream chemistry/geometry gates can still fail and must remain in the denominator.

Raw Ca occupies `B:LIG_B1/CA1`. The desired normalization is residue name only: `B:LA1/CA1`, retaining elemental Ca, atom name CA1 and every coordinate. This matches the existing crystal normalizer's native-element convention. The pinned standard-only protonator expects one residue called LA; the current AF3 normalizer renames only elemental La. Copying the canonical LA1 atom selector would therefore be wrong. The parent is implementing an isolated adapter; this task did not modify a production preparer.

## Source checks

All 250 CIFs, 50 input JSONs, 500 confidence JSONs and 50 ranking CSVs are pinned. Each canonical geometry matches exactly one original La-conditioned sample by atom identities and coordinates. Eight are also byte-identical; the others differ in CIF metadata. All 250 coordinate fingerprints are distinct, with no missing files. Each protein retains its exact sequence and protein atom identities across samples, and ligand inventories match the inputs.

All 25 Ca/La input pairs have identical protein records including MSA/template content, identical PQQ definitions, seed 1, no custom CCD and no bonded-atom pairs. Conditioning changes the metal input and job name. Every output reports the same AF3 build:

`AlphaFold-beta-20231127 (763c1867cead6417bbb16b3bdaf40f3bee822f3997ab248d07626bd1a3b678da)`

These are samples from one AF3 method and seed, not independent folding methods or thermal populations. The historical generic ligand-normalization schema containing “protenix” does not make them Protenix outputs.

## Frozen groups

All five Ca-conditioned samples enter the primary pool. For La, retain the four listed samples and reserve the canonical sample for replay. Selection does not use confidence or energy.

| Reference | Established class | Canonical La sample | Primary La samples |
|---|---|---:|---|
| a0a3f2yly8 | La | 1 | 0, 2, 3, 4 |
| a0acd6b9f2 | La | 3 | 0, 1, 2, 4 |
| c5atj3 | La | 2 | 0, 1, 3, 4 |
| c5axv8 | La | 1 | 0, 2, 3, 4 |
| c5b120 | La | 0 | 1, 2, 3, 4 |
| i0jwn7 | La | 3 | 0, 1, 2, 4 |
| mmol_1770 | La | 3 | 0, 1, 2, 4 |
| mmol_2048 | La | 4 | 0, 1, 2, 3 |
| q88jh0 | La | 0 | 1, 2, 3, 4 |
| q89gy2 | La | 3 | 0, 1, 2, 4 |
| q92wy9 | La | 4 | 0, 1, 2, 3 |
| a8r3s4 | Ca | 0 | 1, 2, 3, 4 |
| atq70401.1 | Ca | 0 | 1, 2, 3, 4 |
| bbl57595.1 | Ca | 2 | 0, 1, 3, 4 |
| o24759 | Ca | 3 | 0, 1, 2, 4 |
| p12293 | Ca | 4 | 0, 1, 2, 3 |
| p15279 | Ca | 2 | 0, 1, 3, 4 |
| p16027 | Ca | 2 | 0, 1, 3, 4 |
| p38539 | Ca | 0 | 1, 2, 3, 4 |
| q4w6g0 | Ca | 2 | 0, 1, 3, 4 |
| q60ar6 | Ca | 1 | 0, 2, 3, 4 |
| q88jh5 | Ca | 2 | 0, 1, 3, 4 |
| q8gr64 | Ca | 2 | 0, 1, 3, 4 |
| q9l935 | Ca | 1 | 0, 2, 3, 4 |
| q9z4j7 | Ca | 0 | 1, 2, 3, 4 |

## Interpretation and reproduction

The 25 proteins and labels are consumed calibration/development data. Additional conformers test structural robustness and conditioning bias; they do not add independent biological validation. The canonical panel already separates by donor composition and charge. Group results by protein and homolog family. Prior energy exposure outside the canonical selection was not globally audited; these are unselected samples, not claimed blind controls.

The earlier 50-source numeric subset proposal in `PQQ_CONFORMER_INPUTS.json` remains an unexecuted planning record. The final source set is **all 250** in `PQQ_ALL250_SOURCES.json`. The parent controls the predeclared aggregation and unchanged scoring recipe. This task launched zero new folds, protonations, energies or jobs.

Rebuild the inventory and source manifest without scientific execution:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/accommodation_controls_20260920/inventory_pqq_folds.py
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/accommodation_controls_20260920/freeze_all_pqq_sources.py
```

This is an archive read/hash operation over about 3.7 GB. No GPU or Slurm submission was used. Baseline and active jobs are unchanged. The parent will prepare and score through a separate finite manifest.
