# Deviations from PLAN.md — two-coefficient predictor test, 9 November 2026

The frozen design was run as written. Additions outside it:

1. A Dy-only RMSE on the 57 complete-series ligands was printed (M1 route 0.650, R2-all 0.744).
2. The R2-all-nou column was added to the per-family table, and the named-ligand residuals were added.
3. One design choice was fixed in PLAN.md before any fit and is the subagent's own: the radius basis
   is normalised to Dy and the penalty grid is 2 × M1's. This makes R2-dy identical to M1.
4. `workspaces/ln_structure_increment_20261030/exclusions.tsv` was left unopened as possibly reserved,
   so the FREEZE pin check covers 11 of 12 files.
5. **Report authorship.** `REPORT.md` and this file were written by the lead agent from the subagent's
   reply and the saved tables, because the harness returns subagent findings as text. The lead agent
   checked the scores against `summary.json` and `per_family.tsv`.
