# Whole-chain OMOL fails canonical calibration

All25 calibration proteins were computed, but the classes do not separate:
min(La)−max(Ca)=**−306.006296kcal/mol**, versus the frozen required>0.02.
No research decision bands were issued. This version is rejected as a general
canonical discriminator; the production ORCA baseline remains unchanged.

| Calibration class | Cases | Minimum R_coord | Maximum R_coord |
|---|---:|---:|---:|
| La-associated | 11 | 35.784842 | 609.423542 |
| Ca-associated | 14 | 29.766735 | 341.791139 |

The La minimum is C5AXV8; the Ca maximum is P38539. Both lie within the dataset's
reported charge range, so the failure cannot be assigned solely to the four
flagged out-of-range proteins. This does not isolate every physical cause.
All whole proteins exceed reported training sizes. No case, charge, geometry,
threshold or water inventory was changed to improve these results.

1H4I and4MAE have raw values101.614126 and110.647960kcal/mol, preserving their
earlier relative ordering. They receive **no calibrated classification**, since
the calibration failed. The raw writer's0/3 transfers reaching a region means
no calibrated decisions were available, not three observed wrong classes.
1KB0 remains unavailable because its whole-chain preparation crosses missing
backbone structure. Its four computed endpoints are retained as invalid
preparation diagnostics. Denominators remain25calibration and3transfers,
with2valid raw transfer scores and1unsupported preparation.

## Numerical implementation and cost

Job1200819 completed104 new energy-only forwards with no execution failure;
eight earlier crystal endpoints were reused after exact input/receipt checks.
Native component accounting passes. Peak memory, each unrounded term, sources,
charge/multiplicity, model/executable identities and all receipts are retained
in the full collection/report.

Actual allocation:3982wall seconds on oneA5000/16CPU,3982GPU-seconds and
63712allocatedcore-seconds; reported CPU4179seconds (scheduler integer-second
precision). Preparation/repreparation and report costs are separate. The final
audited report took345.65wall/242.68CPU seconds,537960KiB peakRSS.
Zero newDFT,solver,training or force calculations. The real integration test
confirms that even successful1KB0 raw endpoints cannot acquire a score/class.

## What this means

Numerical execution is credible under the tested implementation. This version
fails predictive calibration and the independent
[disconnected-spectator consistency test](INTACT_SPECTATOR_REPORT.md).
Its initial five-case improvement remains an observed development result; it
does not transfer to the larger reference set. Whole-chain inference is
computationally affordable at tens of seconds per endpoint, but this scientific
failure prevents promotion.

The [checkpoint audit](CHARGE_EMBEDDING_AUDIT_REPORT.md) identifies global charge
conditioning and exact charge-category degeneracy as concrete representation
limitations. A new common ionic reference cannot fix the raw PQQ reversal,
as the [zero-inference reference screen](SEPARATED_REFERENCE_SCREEN_REPORT.md)
shows. Retain the baseline and pursue a separately declared representation
change; do not recalibrate this failed version into apparent success.

Compact values: INTACT_CANONICAL_RESULT.json. Full audited result:
`workspaces/mace_omol_20260917/intact_panel_report_v1/`.
Replay commands and immutable source versions: INTACT_COMMANDS.md.
