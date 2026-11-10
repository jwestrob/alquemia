# Deviations from PLAN.md — acetate forensics, 9 November 2026

1. **Recursive grep.** The subagent ran a recursive `grep -r` over the repo for about 2 minutes and
   then killed it. It produced no output and left no process. This broke the standing rule against
   recursive traversal on the cluster.
2. **Term split.** `exchange_decomp.py` reads the primary-structure output even when a conformer was
   selected. `legs.py` reads the selected structure's output, so the acetate term split here is new.
3. **NTA value.** NTA is −4.60 here, because it includes Dy_nta_w4. The brief quoted −4.26.
4. **DLPNO scope.** The DLPNO comparison is limited to n = 0 combinations because the DLPNO water job
   had failed. An isomer analysis was added to work around that.
5. **Report authorship.** `REPORT.md` and this file were written by the lead agent from the subagent's
   reply and the saved tables, because the harness returns subagent findings as text. The lead agent
   checked the headline numbers against `legs.json`, `state_treatments.json` and `dlpno_isomers.json`.
