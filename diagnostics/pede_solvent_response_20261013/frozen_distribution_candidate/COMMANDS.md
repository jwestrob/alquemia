# Completed candidate: commands and restart

All paths below are actual Biotite artifacts. The canonical and transfer
calculations are complete; do not resubmit them because an old completion event
arrives or a monitor is gone. Root acknowledged their actual results.

## Read the result

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
cat diagnostics/pede_solvent_response_20261013/frozen_distribution_candidate/REPORT.md
cat workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/summary_v1/SUMMARY.json
```

Joinable source table:
`workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/summary_v1/SOURCES.tsv`.
It includes the exact original/recovery IDs, conditioning metal, expected class,
candidate/baseline contrasts, old-band transfer, selections, works and missing
reasons. Scores use their recorded protocol-specific calibration; no universal
zero or cross-backend aquo offset is applied.

## Reproduce collection without any electronic or component calculation

This optional command uses one Slurm CPU and refuses to overwrite its output:

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
sbatch --nodes=1 --nodelist=node-64-768g-5 \
  --output=workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/recollection_%j.out \
  --error=workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/recollection_%j.err \
  workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/scoring/RUN.sbatch \
  workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/scoring/SCORE.py \
  --config workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/SCORE_CONFIG.json \
  --sha 2f9678683d988a4a983d6f7e2aaabdd312303e3b943cd9cab5bff80a9bcf9985 \
  --stage transfer \
  --native-collection workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/transfer_v1/FINAL_COLLECTION.json \
  --canonical-calibration workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/CANONICAL_COLLECTION.json \
  --output workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/RECOLLECTED_COLLECTION.json
```

This is a reproduction command, not a pending job. An agent submitting it must
arm the standard completion wake on the returned job ID before yielding.

## Original execution and recovery records

Within `workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/`:

- `manifest.json`: pinned native model, inventory, implementation, libraries,
  physical settings and unchanged precision gates; SHA256
  `bcd9fc9d8923139dfcb8d7b8d0191c5c29484e700ddc8ffd488c1dea09c10b41`.
- `CANONICAL_SUBMISSION.json`: actual jobs 4007/4008, exact commands, wake PID3913280.
- `TRANSFER_SUBMISSION.json`: actual jobs 4009/4010, frozen calibration pin,
  exact commands, wake PID3915830.
- `SUMMARY_SUBMISSION.json`: actual job4014 and wake PID3920765.
- `canonical_v1/` and `transfer_v1/`: every native charge replay, cavity, matrix,
  qualification residual and GNU-time receipt; no missing values filled with zero.
- `CANONICAL_COLLECTION.json`, `FINAL_COLLECTION.json`: full scored outputs.
- `CANONICAL_ROOT_ACK.json`, `TRANSFER_ROOT_ACK.json`, `SUMMARY_ROOT_ACK.json`:
  root read/interpretation receipts.
- `FINAL_ACCOUNTING.txt`: actual terminal scheduler records, not estimates.

The independent saved-score regression4006 is under
`workspaces/pede_solvent_response_20261013/frozen_distribution_score_tests_v1/`.
It is complete; its queued notification does not request another test.

## Status and limits

Production is unchanged. This is an archived-pool research scorer, not yet an
end-to-end new-protein command or an automatically selected production mode.
Do not silently use its bands for the released scorer or select its output only
when it agrees with a desired label. A new molecular round needs a declared
physical question within the existing authorization, not a blanket repeat of
these completed manifests. See the final report for the unresolved transfer issue.
