# Re-read the completed two-target sensitivity test

Run from the repository root. The following creates a new aggregate export using
only the genuine completed artifacts; no preparation or molecular execution:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 diagnostics/plm_fold_sensitivity_20260920/collect.py \
 --submissions diagnostics/plm_fold_sensitivity_20260920/SUBMISSIONS.json \
 --prior workspaces/accommodation_torsion_20260920/primary_result_v1.json \
 --scheduler diagnostics/plm_fold_sensitivity_20260920/SACCT.txt \
 --output workspaces/plm_fold_sensitivity_20260920/replay_v1.json
```

Original full result: `workspaces/plm_fold_sensitivity_20260920/RESULT_v1.json`.
Each standard plan/request/result/ensemble and source preparation is retained under
its named target directory. `SUBMISSIONS.json` contains exact submitted commands
and wrapper hashes. Do not resubmit completed calculations into these immutable
paths. The wrapper preserves failed preparation as an unavailable group and never
selects a replacement fold.

For a separately authorized new request use the existing
[standard three-fold instructions](../pqq_ensemble_20260920/COMMANDS.md). The
isolated PLM test does not authorize a production scan or change defaults.
