# Deviations from PLAN.md — four-carboxylate rule on predicted structures, 9 November 2026

The frozen design (PLAN.md sha256 `e16399ac…f4b2`) ran as written. Notes:

1. **Compact EF-hand group.** Added at the lead agent's request after the freeze; descriptive and
   outside the frozen reading.
2. **Partner tier.** Taken from the benchmark's `site_ranks.tsv` by exact site-id match (256 sites in
   28 proteins); it was not re-scored here.
3. **LanD site.** The 9C90 holmium site used for the reference table is symmetry-bridging.
4. **Reference sequences.** Precursor sequences were folded; the reference table gives no mature
   sequence.
5. **Transfers.** No site transfer failed; 16 had a local RMSD above 2 Å.
6. **Tracked files.** Candidate identifiers and sequences come from unpublished datasets. They are
   left untracked (`.gitignore` in this directory), and `REPORT.md` gives group-level numbers only.
7. **Report authorship.** `REPORT.md` and this file were written by the lead agent, because the
   harness returns subagent findings as text.
