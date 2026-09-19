# Completion verification

- All 25 cases completed for DFT and each MACE implementation. The finalizer
  ran successfully against the actual artifacts, with 50 DFT and 100 MACE
  endpoint receipts, four failed MPI startups, and all terminal job records.
- Both MACE implementations give 25/25 correct calls and exactly reproduce
  all frozen scores. No fitted bands, physical input changes or new model.
- Archived 25+3 reference replay retains 1KB0 in the unavailable denominator.
  Two supported crystals pass; calibration and consumed transfer stay distinct.
- Prepared-input-to-score times include validation, loading and reporting.
  Shared DFT results are counted once in the combined costs. Raw preparation,
  queue waiting and unknown local development totals are not misreported.
- The planned speed gate passes; the combined qualification stays false.
  Q9L935's DFT score exceeds reproduction tolerance while retaining its class.
  The separate La boundary crossing stays visible. Neither issue is concealed.
- Real-fixture source-selector, prepared-interface and XYZ-cache tests each
  passed four tests (14.902, 28.420 and 18.552 seconds respectively). Their
  existing logs were inspected; no fabricated scientific fixtures were used.
  The earlier two archived-score/parser tests also passed. Reporting was
  exercised both against running jobs (correct refusal) and completed jobs.
- Baseline/default unchanged; no new scientific protocol, score threshold,
  production rescore or push. Latest accuracy-first steering is recorded.

Machine-readable completion audit:
`workspaces/mace_pqq_utility_20260918/completion_audit_v1.json`.
Finalizer execution receipt: `finalize_execution_v1.json` in the same workspace.
The scientific outcome and limitations are in [REPORT.md](REPORT.md); commands
for replay without model calls are in [COMMANDS.md](COMMANDS.md).
