# Four-parameter radial continuation: independent review

**No release blocker found.** This review covers the changed start, additional budget, immutable native-cache admission and reduced full-result checkpoint frequency. It does not qualify the force field physically or release donor transfer.

## Exact continuation and reuse

The starting vector, eight native energies/forces, objective0.4970127555469721 and original-center diagnostics replay the selected17589 result exactly. Original bounds, fixed epsilons, training data, body definitions, objective and direction helper remain unchanged. The17576 data baseline and original SI bound origin remain distinct from the new search start.

The new ledger allows32 additional outer steps,384 new vectors and4096 new program invocations. Parent actual costs remain12 steps,132 vectors and1344 programs. The additional scope does not reset or rewrite those records.

All672 complete parent native cells are indexed by effective physical input. Eight starting cells undergo exact raw E/F replay without native calls. Inherited results retain their original native authority, implementation and cache key; the current request identity is separate. Exact opposite-metal reuse passes the author's narrow controls. Baseline17578 provenance correctly leaves `native_implementation` null while retaining its actual native authority. Protected parent result/config/ledger pins remain unchanged.

## Checkpoint and failure behavior

The actual nested `batch` function was extracted from frozen REFINE.py and exercised with two fake workers, without engine imports or native calls:

- Successful batches save only at start and end; duplicate vectors reuse one entry.
- The start snapshot declares every scheduled point and all eight cells before execution.
- A deliberately failing callback saves immediately while its sibling is still running.
- The sibling's completed outcome is retained after the failure. The final failed point preserves its partial cell progress and the full denominator.
- A failed batch raises before selection or another batch; the global stop is set. Native per-cell results and evaluator ledger persistence remain unchanged.

The success simulation produced2 full writes; the failure path produced4: start, immediate failure, drained batch, detailed failure. This tests evidence retention under handled worker failures; it is not a fault-tolerant restart claim for process or node loss. Selection still uses complete batches, unchanged ordering and current-wins ties.

## Evidence and cost

Final independent result: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_radial_continuation_independent_review_v1/REVIEW_v2.json`, SHA `cad00686fa4d15f1229c252fa0147ff189337ebcf85d2b691992216a10dfd527`.

71 saved checks passed in1.028 seconds with one numerical thread and at most2 fake callback workers. Zero native calls or reserved-validation outcomes. The author's32 cache/budget controls are reused rather than repeated. Initial review receipt/source remain preserved; v2 fixes only the review's check-count bookkeeping by keeping the fresh-output guard outside the scientific check list. Total local wall time for both review executions was2.061 seconds.

| Evidence | SHA256 |
|---|---|
| REFINE.py | `f6550a635a6f4fd0ee3e801011bb99c869d4ca7352f209c656fee52d59ca96cd` |
| EVALUATE.py | `86ad2f4c36afc441d0338a92a4aa3eb5ca6d1477ed2d7e57d1a11bcbbced9fce` |
| SPEC.md | `95e29d64fa557669554d39d9f4a07abe89aaf2eee96ac030b0d414a6405ab9bf` |
| Evaluator config | `c6aa9d2ad27a60027364e7c86cc9b68ac62475d58281cc3c777977864ff0583e` |
| Actual preparation_v2 | `523bc150b2d92665ff3d16c9b9881a9544947f7140c7a50c5a4ea98a420e1413` |
| Author saved_controls_v2 | `e3464b87df385479f7f0dcde5f27b02a9aa51ef7a983559d8d77bbe268eabd06` |
| Reviewer source | `2d8780c29127f1c6c05329df1c85e43eeecf0ceabc05c754eee0f9b42d6f697e` |

Root retains submission authority. Passing development controls would support representability on these consumed fixtures, not bulk-water behavior, protein transfer, affinity accuracy or electronic-state qualification.
