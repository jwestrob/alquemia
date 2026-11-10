# Deviations from PLAN.md — lattice-strain test, 9 November 2026

1. **Feature source.** `workspaces/ln_selectivity_calibration_20261027/features.tsv` has no charge or
   denticity columns. `q_lig` and `nlig_La_cn8` were taken from `surrogate_nist.tsv`, which derives
   them from MACE-built complexes.
2. **Reading of F1.** F1 was implemented as a two-parameter per-ligand fit on the leave-one-ligand-out
   SVD shapes. The existing regression map from `series_rank.py` (39 coefficients shared across
   ligands) is reported beside it as F1r. The "equivalent to SVD" reading passes against F1 and fails
   against F1r on held-out metals (+0.043).
3. **Added after the first run.** The projection of the SVD shapes onto the Δr, Δr² plane, and the
   cubic and Gd-step breakdown of the per-metal residual, were added once the scores were known. They
   are descriptive and carry no frozen reading.
4. **Crash and rerun.** The first run crashed at the feature merge after the scores had printed. The
   rerun reproduced the same scores.
5. **Report authorship.** `REPORT.md` and this file were written by the lead agent from `run.log` and
   `summary.json`, because the analysis subagent's tools declined report-file writes. The lead agent
   checked the summary numbers against the subagent's reply before writing.
