# Nikasha: bounded structural increment, 2026-10-30

Completed. **Target missed; retain donor counts plus charge and end automatic expansion.**

- [Results and interpretation](REPORT.md), [machine-readable results](RESULTS.json)
- [Frozen plan](PLAN.md), [groups](groups.tsv), [splits](splits.json), [pins](FREEZE.json)
- [Disclosed repairs and failed attempts](DEVIATIONS.md), [actual costs](COMPUTE_COST.json)
- [Per-family ablations](per_family.tsv), [paired improvements](paired_improvements.tsv)
- [Published reference per family](published_per_family.tsv), [descriptor sensitivity](descriptor_sensitivity.tsv)

Row-level NIST-derived records, exclusions, predictions, complete overlap inventory,
coordinates and all attempts are under `workspaces/ln_structure_increment_20261030/`.
They are deliberately not redistributed in this public repository, following the
existing NIST pipeline convention. The report links each table directly.

From the repository root on biotite:

```bash
bash diagnostics/ln_structure_increment_20261030/reproduce.sh
```

This uses a new output directory and existing pinned assets/environments. It
submits only the declared statistical/MMFF/graph-inference stages with
`--qos=standard_expand`, then waits for saved summaries. The original NIST source
is `workspaces/ln_selectivity_calibration_20261027/nist46_dat/`; existing extraction
and historical donor SMARTS are reused. See the original calibration README if
that local source must be recovered from the NIST installer.

The compatible published-model environment is pinned in `chemprop-lock.txt`.
The analysis environment is recorded in `analysis-environment.json` and uses the
existing `/groups/banfield/users/jwestrob/envs/lnsel/bin/python`. Asset URLs and
Git blob hashes are in `assets.json` / `additional_assets.json`; local SHA256
receipts are in the workspace's upstream manifests and `LOGKPREDICT_ASSETS.json`.

No scientific choices should be changed in a reproduction. The grouping repair
is post-initial-results and explicitly flagged; the legacy LOGKPREDICT numeric
example remains failed. Overlapping pretrained scores are not independent tests.
