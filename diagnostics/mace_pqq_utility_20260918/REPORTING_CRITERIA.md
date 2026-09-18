# Reporting conventions before the complete timing panel

GOAL.md specifies at least1.5x faster median prepared-input-to-score time and
lower total case time. The acceptance metric is therefore
`median(DFT case seconds) / median(MACE case seconds)`, with all25cases included.
The median of the25matched pair ratios is also useful and is reported separately;
it is not substituted for the originally stated median-time criterion.
Early progress updates quoted paired medians and were explicitly provisional.
Whole-allocation elapsed-time ratio is included to expose startup overhead.

This clarification was recorded while only6DFT/10MACE cases were complete. No
new calculation, case selection, scientific parameter, calibration or threshold
was introduced. The frozen0.01score-reproduction tolerance and literal historical
classifications remain separate. Report a failure or a numerical boundary crossing
explicitly; do not adjust a biological decision band to force a passing result.
