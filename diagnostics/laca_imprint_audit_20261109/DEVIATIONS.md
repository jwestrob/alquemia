# Deviations from PLAN.md — La/Ca imprint audit, 9 November 2026

1. **Environment.** Geometry ran under `conda_envs/lanm_qmmm` (gemmi 0.7.5) because `envs/lnsel` has
   no gemmi. HTTP work and statistics used `lnsel`.
2. **Deposit count.** PLAN.md says the 125 Ca sites come from 42 deposits. The count is 34; 42 belongs
   to the full 175-site panel. The power calculation uses 34.
3. **Chain copies.** Within one entry, chain copies are grouped by their set of coordinating residue
   numbers. The bootstrap unit remains the protein.
4. **Engineering additions.** An alignment cache and a guard for residues with no `label_seq` were
   added. Neither affects outcomes.
5. **Heavy / light split.** The split by lanthanide was added after the primary result and is
   descriptive; it carries no frozen reading.
6. **Report authorship.** `REPORT.md` and this file were written by the lead agent from the analysis
   subagent's reply and the saved tables, because the harness returns subagent findings as text. The
   lead agent checked the primary statistics against `results/stats_summary.json` and
   `results/per_protein.tsv`.
