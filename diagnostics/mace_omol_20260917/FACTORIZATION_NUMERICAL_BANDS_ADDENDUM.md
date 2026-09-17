# Execution-specific numeric calibration records

Declared2026-09-17 while canonical job1200830 runs; its final calibration has
not been inspected or issued. No change to that job or its frozen decision rule.

The two-call and four-call formulae agree to roundoff, but inclusive extrema
used as bands can differ by tiny floating-point amounts. Do not silently apply
one execution's exact boundary to another score and call a roundoff change a
biological disagreement. Do not tune an epsilon after seeing a borderline case.

If the full-panel factorization check passes the predeclared0.01 tolerance,
produce a separate, explicitly typed two-call calibration record using the same
25 designated calibration cases, their actual saved bound endpoint values and
the fixed1H4I disconnected model terms. Compute the same predeclared extrema
and require the same gap>0.02. No new labels, subset selection, fitted coefficient
or reference-anchor choice. Transfer cases never set the boundaries.

Preserve the original four-call report/bands. Record a score_evaluation ID and
require it to match the supplied calibration. Two-call scores use their matching
numeric reference; four-call scores use theirs. Both represent the same tested
masked descriptor within the established numerical tolerance, not different
scientific hypotheses. Failure of equivalence blocks the two-call classifier.
No non-PQQ affinity band, physical zero, baseline refit or default promotion.
