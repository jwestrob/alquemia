# Grouping terminology correction — 2026-09-23

The hash-pinned REPORT.md mistakenly calls the La4, Ca5 and three-of-four
summaries arithmetic means. The executed implementation actually uses the
predeclared **strict median** of the available member scores, requiring every
specified member. The balanced summary is the **equal arithmetic mean of the
complete La4 and Ca5 medians**.

The code, actual numerical scores, decision counts, reference and missing-member
policy are unchanged. La4 remains23 correct/2 unavailable; Ca5 remains22/3;
balanced21/4; and the overlapping La triples94/6. All available grouped calls
remain correct. These are correlated robustness summaries, not independent
biological observations. No scientific execution or recalibration was performed.

The exact path is `standalone_xtb_transfer.compare -> nikasha_pool_compare.aggregate
-> accommodation_fold_proposals.strict_summary -> statistics.median`.
The balanced-arm arithmetic is explicit in `standalone_xtb_transfer.compare`.
The original report is retained byte-for-byte so its artifact pins remain valid.
