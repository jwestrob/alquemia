# Motion-envelope pilot: prepared inputs and actual execution

Run from the repository root. `CPU` below is the existing chemistry environment;
all paths are arguments. The source selection and preparation are immutable.

```bash
CPU=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
$CPU scripts/motion_envelope_run.py validate \
  --manifest workspaces/motion_envelope_20260923/pilot34_searches_v1/manifest.json
$CPU scripts/motion_envelope_scalar.py validate \
  --manifest workspaces/motion_envelope_20260923/pilot34_scalar_origins_v1/manifest.json
$CPU -m unittest discover -s tests -p 'test_motion_envelope*.py' -v
```

The saved `SUBMISSION.json` files contain the exact successful scheduler commands,
manifest/wrapper hashes and job IDs. Origin MACE1211294 completed68/68. Search
1211342 completed all 68 starts; strict-origin1211343 completed 136 scalar cells.
These are finite tasks: inspect their terminal collections before doing anything
that would submit them again. Execution uses each stage's frozen implementation.
The wrapper collects actual terminal results even when the executor fails.

Final assembly joins actual strict scalar origins to the separately pinned
native-only search manifest. It does not rewrite search inputs or insert missing
components. The prepared common pool holds exactly origin/Ca-proposal/La-proposal,
with numerical coordinate deduplication and both metals scored on every member.

```bash
$CPU scripts/motion_envelope_pool.py prepare_pool \
  --proposals workspaces/motion_envelope_20260923/pilot34_searches_v1/collection.json \
  --scalar-origins workspaces/motion_envelope_20260923/pilot34_scalar_origins_v1/COLLECTION.json \
  --agreement diagnostics/motion_envelope_20260923/EXECUTION_PLAN.md \
  --output workspaces/motion_envelope_20260923/pilot34_pool_v1
```

The output directory must be new. Pool preparation emits the finite candidate
scalar manifest using the exact qualified native TolE1e−10 fresh profile. No loose
scalar cache is accepted. Actual cross-MACE and scalar execution remains allocated
work, using the existing runners. The final report records subsequent job IDs.

Only the25 designated canonical rows calibrate the new reference; three crystals
and six probes are held out. The matched strict-tenfold baseline uses
`strict_native_pool_20260923/run_v1/COMPARISON.json` plus
`strict_native_comparator_20260923/run_v1/COLLECTION.json`. Previous loose scores
remain historical columns. Full100-triple molecular transfer is not part of this
pilot's execution scope.

Cross-MACE1211421 completed68 cells. Candidate scalar1211422 completed 272 cells
with the same strict profile. Collect the terminal pool and form the two separate
canonical-only references (adaptive variants and the authorized origin ablation):

```bash
$CPU scripts/motion_envelope_pool.py collect \
  --manifest workspaces/motion_envelope_20260923/pilot34_pool_v1/manifest.json \
  --output workspaces/motion_envelope_20260923/pilot34_pool_v1/COLLECTION_FINAL.json
$CPU scripts/motion_envelope_compare.py reference \
  --collection workspaces/motion_envelope_20260923/pilot34_pool_v1/COLLECTION_FINAL.json \
  --output workspaces/motion_envelope_20260923/REFERENCE_v1.json
$CPU scripts/motion_envelope_compare.py compare \
  --reference workspaces/motion_envelope_20260923/REFERENCE_PINNED_v1.json \
  --strict-primary workspaces/strict_native_pool_20260923/run_v1/COMPARISON.json \
  --strict-supplement workspaces/strict_native_comparator_20260923/run_v1/COLLECTION.json \
  --strict-reference workspaces/strict_native_pool_20260923/run_v1/REFERENCES.json \
  --strict-origin-reference workspaces/strict_static_ablation_20260923/run_v1/REFERENCE.json \
  --triple-reference workspaces/union_triple_pilot_20260923/REFERENCE_v1.json \
  --triple-transfer workspaces/union_triple_transfer_20260923/COMPARISON_v1.json \
  --output workspaces/motion_envelope_20260923/COMPARISON_v1.json
$CPU diagnostics/motion_envelope_20260923/summarize_costs.py \
  diagnostics/motion_envelope_20260923
```

These exporters create new immutable JSON files; use the saved completed files for
review, or a new output name for an independent replay. Never overwrite an existing
scientific result or resubmit a completed molecular stage to refresh a report.

The completed comparison uses `REFERENCE_PINNED_v1.json` (SHA256
`4a1d3ae6e83b5b6ace408f729c39390d44f673554e338df932fd38fc971f67d9`).
It preserves all scientific fields and the frozen timestamp of `REFERENCE_v1`;
its implementation pin resolves to the byte-identical archived reference writer
in `analysis_reference_v1/`. The original export is preserved. A later component-key
alias repair affects reporting only, not calibration or molecular results.

Final compact outputs are `RESULT.json`, `PROBES.csv`, `REPORT.md` and
`COSTS_v1.json` in this directory. `export_report.py` records their exact source
artifacts; its already completed invocation is retained in `EXPORT_v1.log`.
Do not regenerate over these outputs. The five completed actual-fixture suites
record 18 passing tests in `TESTS_{preflight_v2,execution_v1,join_v1,pool_v1,final_v1}.txt`.
