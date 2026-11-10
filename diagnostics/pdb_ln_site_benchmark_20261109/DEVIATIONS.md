# Deviations from PLAN.md — PDB lanthanide-site benchmark, 9 November 2026

The frozen design (PLAN.md sha256 `c9dc815c…7853`) ran as written. Notes:

1. **Donor-origin labels.** The lead agent's request to classify each site as single-chain,
   inter-chain or symmetry-bridging arrived before the freeze and is part of PLAN.md. The frozen
   primary reading uses all sites; the single-chain version is co-reported.
2. **Choices fixed in the plan that shape the result.** Each chain is scored in isolation. Partner
   deposits extend to 3.0 Å resolution. Copies within the biological assembly are unresolved, so every
   donor from another crystal copy counts as symmetry-bridging.
3. **Dropped material.** 33 of 879 labelled sites fell out of the stripped tier through the
   four-chain cap or chains under 20 residues. 223 partner transfers failed the eight-Cα-pair rule.
   The existing scanner (B1) errored on 8 of 1,286 scored units.
4. **Code check before scoring.** One smoke run on two structures outside the benchmark (4ICB, 1EXR);
   nothing was tuned.
5. **Development split.** 40% of the 30%-identity clusters were set aside for development and never
   used, because v0 has no tuned parameter; results are reported for all clusters and for the
   held-out test split.
6. **Report authorship.** `REPORT.md` and this file were written by the lead agent, because the
   harness returns subagent findings as text. The lead agent read `results/stats.log` for every number.
