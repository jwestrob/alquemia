# PedE diagnosis complete — 13 October 2026

Read [the integrated report](REPORT.md). Jacob asked for an investigation and
opinion before more experiments. This round used saved artifacts only; no new
molecular energies, folds, optimization, fit, promotion or production edits.

All five rebuilt two-Ca PedE sources and five prior raw La sources were inspected;
four prior La pools are method-compatible, while old La1 remains unscored.
Fine donor geometry, cofactor pose, confidence and proposal motions do not expose
a shared wrong-only defect. The correct3→wrong1 origin-score difference 21.1062
kcal/mol decomposes into −0.1504 MACE and +21.2566 solvent. Saved solvent charges
show much greater PQQ/context redistribution in the wrong sources. All recorded
SCF/component gates pass; physical correctness or a discrete branch switch is
not established. Shared fine cofactor geometry is not experimentally qualified.

The two wrong origins were already inconclusive. Accommodation narrows the
canonical gap 10.2977→3.4722 kcal/mol; wrong1/2 lose only 1.230/2.014 kcal/mol of
Ca-band margin after accounting for their own calibrated scale. The full125
outcomes stay static117 correct/1 wrong/7 inconclusive versus accommodated119/3/3.
No cases, labels or thresholds were changed. These are consumed correlated
structural samples, not fresh blind biological observations.

Jobs3338/3339/3340/3341/3344/3345/3348 are terminal; original reader/export
failures are preserved with explicit successors. Total12 allocated CPU-seconds
at scheduler integer resolution, zero GPU/new molecular calls. Every completion
has a root result-read acknowledgement. Delayed wakes require acknowledgement,
not resubmission. Agents finished their scoped reports and vault notes.

Recommendation for discussion: an independent matched solvent/charge-response
comparison on saved correct3/wrong1 origins, both metals and both media, before
more folding/search. This is proposed only; no automatic molecular follow-on
is queued or authorized by a diagnostic completion event. Preserve existing
production and unrelated scheduler-listed jobs; readiness probes remain off.

Next action: discuss the report with Jacob. Reproduction commands and exact
inputs are pinned in the workspace SUBMISSION/config/manifest records linked
by EXECUTION.json. Do not reopen closed LanM or PQQ campaigns from older notes.
