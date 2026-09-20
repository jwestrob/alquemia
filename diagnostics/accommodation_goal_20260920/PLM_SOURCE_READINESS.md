# Two existing PLM development cases: three-fold inputs located

The two already-diagnosed compressed sites have all three saved La-conditioned
AF3 models. Their exact source sequences, residue numbering, prior homologous
core roles, protein/La/PQQ assembly and fold-input identity are compatible with
declaring a three-fold source request. No new folding is needed.

This is input readiness, not successful chemical preparation or prediction.
No new protonation, carve, energy calculation or PLM rescore ran here. Neither
target has a biological La/Ca label. Both were selected for the existing donor
compression investigation, so they cannot become fresh validation cases.

## Fixed input scope

- `PQQSEQ_07ab500e3df76b30d71c`: all saved seed101 samples0/1/2.
- `PQQSEQ_83440678cbbd658047c9`: all saved seed101 samples0/1/2.

The prior selected-sample record supplies the five homologous residue identities.
Each identity was checked against every actual source model and the original
AF3 protein sequence. The source La atom and PQQ residue were checked directly.
No coordinates or atoms changed. Full source/input/selection hashes are retained.

The first target's sample2 has `admission_pass=false` in the previous single-fold
selection record. It remains included and flagged. This request does not override
that historical decision or assert that the new preparer will support it. A
missing/unsupported member must leave the declared three-fold descriptor unavailable;
do not substitute another sample or report a selected two-fold median.

## Artifacts and reproducible read-only operation

`workspaces/accommodation_goal_20260920/plm_source_requests_v1/` contains two
source-request JSON files and `inventory.json`. The requests use the released
source configuration and developmental ensemble schema. They contain no assigned
expected class and no calculated score.

The actual operation, reproducible into a fresh output directory, is:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 diagnostics/accommodation_goal_20260920/prepare_plm_sources.py \
 --candidate-manifest workspaces/plm_xoxf_all_fixed_core_20260916/retry_1200785/prepared_batch/candidate_manifest.json \
 --release-sources diagnostics/pqq_fast_release_20260920/SOURCES.json \
 --targets PQQSEQ_07ab500e3df76b30d71c PQQSEQ_83440678cbbd658047c9 \
 --output workspaces/accommodation_goal_20260920/plm_source_requests_replay
```

Actual result: both groups and all six raw models verified. Chemical preparation
and energy execution remain `not_run`. The native DFT donor-response experiment
is separate and uses the original selected source, not these new request files.
