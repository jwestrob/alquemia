# Both interrupted MxaF sources recovered

Jobs **2869** (native proposals) and **2873** (corrected solvent components)
completed successfully. Both P12293 Ca-conditioned sources, samples 1 and 4,
now have complete shared geometry pools and classify **Ca**, correctly, under
operational, mathematical-minimum and origin policies with the frozen bands.

This restores coverage lost to the earlier node failure. It does not change the
model or establish an additional accuracy improvement. The original 125-source
result remains **105 correct / 1 wrong / 1 inconclusive / 18 unavailable**.
The explicitly separate recovery overlay is **107 / 1 / 1 / 16**. The other
sixteen unavailable sources and both non-correct available outcomes remain visible.

## What ran and what changed

Both La optimizations restarted at the original q0 under the unchanged four-mode
SLSQP policy. They completed in 16 and 14 iterations. Existing Ca proposals,
origins, eight native cells and sixteen solvent inputs were reused. The final
operational choices retain each metal's own proposal for both sources; common-pool
cross-scoring was complete. Both metals lower their local composite energy, while
both paired scores remain on the Ca side of the unchanged calibration.

Actual Slurm checks passed:

- Exact task/selector reconstruction and all four q0 physical coordinate maps.
- Both original CPU runtime controls: zero reported energy error and maximum
  force discrepancies below 1e-13 eV/Å.
- Both new admitted La proposals and opposite-metal cells.
- Exact replay of eight old solvent components and both origin scores before
  new solvent calls.
- All required corrected xTB/tblite component, population, SCF/state and charge
  checks; no unreceipted interrupted new calls.

No new Ca optimization, origin evaluation, DFT, folding, calibration, chemistry
change or production promotion occurred. The failed/interrupted primary records
and their costs remain preserved. All scientific preparation/checks/evaluations
ran inside Slurm; reporting only reads their completed artifacts.

## Measured incremental cost

| Stage | Actual work | Slurm elapsed / CPUs | Process measurements |
|---|---|---|---|
| Native recovery | 35 actual forwards, including two controls and two cross cells | 441 s / 7 | 2330.266 CPU-s across calls; warm worker peak RSS 7,208,096 KiB |
| Solvent completion | 8 primary xTB and 8 independent tblite calls | 105 s / 12 | 397.298 CPU-s; maximum backend RSS 109,452 KiB |

Application receipts also retain preparation/collection timing and allocated
CPU-seconds separately from process CPU time. These are measured costs for these
two sources, not a promised throughput or resource-utilization measurement.

## Durable outputs

- Native: `workspaces/ca_coverage_recovery_20261012/mxaf_v2/FINAL_COLLECTION.json`
- Solvent: `workspaces/ca_coverage_recovery_20261012/mxaf_solvent_v1/FINAL_COLLECTION.json`
- All125 primary-preserving comparison: same solvent directory,
  `PRIMARY_RECOVERY_OVERLAY.json`.
- Exact old-input replay: same directory, `REUSE_CHECK.json`.
- Compact source/result/cost pins: `RESULT.json` beside this report.

This finite technical-recovery branch is complete. Do not resubmit either job.
The separate two-Ca refold/spectator-preparation hypothesis addresses a different
subset and cannot retroactively replace these results.
