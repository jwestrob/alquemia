# PLM all-XoxF expansion: active launch

On 17 September 2026, Jacob’s approved all-XoxF AF3/current fixed-core workflow was resubmitted after a host-launcher compatibility repair.

- MSA **1200793** is actively searching UniRef30 on `node-64-768g-15`, using all 64 allocated CPUs. Native ColabFold preflight passed on this node with Apptainer 1.1.3; actual search and ungapped-prefilter progress are recorded.
- AF3 **1200794** → preparation **1200795** → La/Ca endpoints **1200796** are queued with immediate-predecessor `afterany` dependencies. Valid available subsets proceed; every absent or failed target is explicitly retained.
- All **176** exact full proteins remain selected; two completed current-method pairs are reused. There are 169 new MSA queries, 174 AF3 monomers × 3 samples, and up to 174 new La/Ca pairs. The frozen fixed-core chemistry, selection, gates and calibration are unchanged.
- Watcher PID **1692026** monitors all four jobs and queues notices to the original PLM session; its first start notice was successfully queued.

Current PLM attempt:
`/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/xoxf_all/retry_1200785/`

Read its `README.md`, `job_status.json`, `launch_receipt.json`, `launch_validation.json` and eventual `completion.json`/`results.tsv`. Input provenance is in the parent `inventory/` directory. Preparation and endpoint products are under `workspaces/plm_xoxf_all_fixed_core_20260916/retry_1200785/` in this repository. This note accompanies the committed preparation adapter and its 16 passing checks; folding had 14 checks plus independent equivalence review.

Initial jobs 1200785/1200786 failed before search/inference because a login-host container-launcher hash was required on compute nodes. The scientific containers/weights were intact. Jobs 1200787/1200788 produced accounting only: two prior scores and 174 missing upstream folds. Those 174 historical rows were mislabeled unsupported; the parent `initial_attempt_review.json` corrects them to workflow failures. All raw files remain preserved. The retry separates host-launcher compatibility from scientific payload identity and distinguishes workflow failure from chemical exclusion.

Current run-plan SHA256: `b05838075b10e15566649da4aaa8a105e1ba7194d2aba72b7bef94bb1cdeedba`.
Folding-manifest SHA256: `cd4d193e10834832c95be030c1f6f118b0cbf4fdc7cf6820aa478280f6147ad7`.

This does not restart cancelled old queue1196402, change the production method, rerun references or import the parallel MACE pilots. The two reused scores remain AG41 17.21776065297137 (indeterminate) and Rokubacteriales 29.971263921125484 (La-supported) kcal/mol. Predicted-structure scores do not establish physiological metal use or measured affinity.
