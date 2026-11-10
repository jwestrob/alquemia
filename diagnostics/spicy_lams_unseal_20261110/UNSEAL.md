# SpyCI-LAMBS outcomes unsealed for 600 of 616 orthologs — 10 November 2026

**Authorized by Jacob on 10 November 2026:** "Unseal the data as you like; just mark that you did it."
Done by the lead agent (Claude Fable) at 2026-11-10 08:38 UTC with `extract_nonreserved.py` in this folder.

## What was opened

Workbook `Jacob/spicy_lams/41589_2026_2176_MOESM3_ESM.xlsx` (sha256 `39193323…99f1`), sheet `Data`:

- raw desorption by ICP-MS, `T:BL` (15 elements × 3 replicates), and
- normalized logD mean and SD, `BT:CW`,

for the **600 rows outside the reserved panel**. The extract is
`workspaces/spicy_lams_unseal_20261110/nonreserved_profiles.tsv` (hash in `EXTRACTION.json`).

## What stays sealed

1. **Every outcome cell of the 16 reserved proteins** of
   `spicy_lams_8family_double_difference_sentinel_v1` (workbook rows 17, 29, 128, 130, 173, 201, 210, 373, 374,
   383, 384, 399, 400, 417, 515, 519). The extraction skips those rows by row index; only their identifiers,
   sequences, taxonomy and total-metal QC columns (`A:S`, `BM:BS`) were written, to
   `reserved_metadata_only.tsv`.
2. **The published cluster and UMAP columns (`CX:CZ`) for all rows.** The extraction never requests them.
3. The article's results sections and the derived preprint stay unread by this project's agents, because they
   name individual orthologs.

## Consequences

- Any statement that no SpyCI-LAMBS outcome has been opened is true only up to 10 November 2026.
- **Sentinel v1** (`spicy_lams_matched_homolog/PREREGISTRATION.md`) is unchanged: its fixed-cluster energetic
  predictions were never locked and its 16 rows are still sealed. It remains evaluable only by a method locked
  without any use of the 600 opened rows.
- A sequence-based protocol (v2) will be hash-locked from the 600 opened rows before the 16 rows are read, and
  scored once with the v1 statistics.
- Analysis plan frozen before any outcome value was read:
  `~/Documents/Berkeley/Lanthanides/lanm_ortholog_profiles_2026-11-10/PLAN.md` on Jacob's Mac
  (sha256 `155dccc5…0363`), copied beside this file.
