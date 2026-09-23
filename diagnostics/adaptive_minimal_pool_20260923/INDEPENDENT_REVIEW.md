# Independent review: minimal adaptive replay and two-source recovery

Reviewer: water_basins, 2026-09-23. Read-only inspection of existing code, physical
manifests, actual completed outputs and comparison algebra. No molecular calls,
new scores, threshold fitting, source changes or shared-code edits.

**No blocking scientific or implementation issue found in these actual results.**
The reduced pool is a justified computational simplification on consumed data;
the two recovered rows improve coverage, not two wrong static predictions.

## Verified

- Both metals uniformly receive exactly origin, adaptive_Ca and adaptive_La.
  Complete six-cell matrices are required. Mathematical minima remain distinct
  from the unchanged 0.1 kcal/mol origin-retention rule. Labels do not enter
  geometry selection. Dropping older candidates changes eleven scores but no
  available class decision in the archived225 replay.
- The separate reference uses only the original25 designated canonical members.
  Its class extrema exactly reproduce their actual frozen rows; crystal/PLM and
  noncanonical folds do not enter calibration. The original and recovered
  comparisons pin the same reference. No threshold was refitted after recovery.
- The actual recovery run uses the unchanged paired q0 force selector, same
  four angular coordinates and scaled SLSQP implementation/settings. Exact
  origin components, source nuclear pairs, charges, physical mappings, model
  and native recipe are reused. Both completed candidates per source are
  required, and each cross-cell changes only the metal/state as declared.
- All16 newly executed native outputs reproduce the saved energies exactly;
  actual receipts confirm task/manifest identity, return code, convergence and
  normal termination. Charge audits pass and paired vacuum/ALPB parameter
  exports match. The MaxIter500 qualification is present. These checks do not
  import the separate continuation experiment or establish unique SCF branches.
- All four final geometries meet the declared physical checks. Maximum heavy
  displacements are P12293 Ca0.3321/La0.5628 Å and Q60AR6 Ca0.8000/La0.5211 Å.
  The Q60AR6 Ca proposal is boundary-constrained, not an unconstrained minimum.
- Actual recovered225 has225 unique rows,75 strict source groups and100 fixed
  three-fold groups. All earlier method fields remain exactly unchanged. Only
  the two specified previously unavailable rows receive a separate new result.
  Missing members still invalidate strict group summaries. No empty group,
  favorable subset or origin-only fallback becomes a successful adaptive score.

## Coverage and interpretation

| Comparison | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
|Released static, all225|203|2|2|18|
|Minimal adaptive replay, all225|202|1|1|21|
|Minimal plus completed recovery, all225|204|1|1|19|
|Released static, matched206|202|2|2|0|
|Recovered adaptive, matched206|204|1|1|0|

The newly available P12293 Ca-sample4 and Q60AR6 Ca-sample1 both remain correctly
Ca-like, as their static origins already were. Historical unavailable records
remain nested beside their new results. The distinct Q60AR6 La-sample0 optimizer
failure remains unavailable; no retry or disguised fallback occurred. The
17 preparation exclusions and one missing native origin also remain in the
225 denominator. These structural repeats are consumed development data from25
reference groups, not225 independent biological labels. Available aggregate
accuracy was already correct; the improvement is individual-fold robustness.

The40% nominal reduction (12 versus20 solvent cells for three versus five
geometries) is arithmetic, not a measured40% production speedup. Reused
optimization/historical costs remain real. The native-solver and boundary
limitations remain material. This review supports the stated result and
coverage accounting, not automatic default promotion or affinity interpretation.

## Exact reviewed artifacts

- `scripts/adaptive_minimal_pool.py` — SHA256 `fdf1b4856afc9abea1925b4781161094ae3dba84264fc85a497ae6d4ac35567c`.
- `scripts/adaptive_origin_recovery.py` — SHA256 `3ae091652f2b358eb15e773d33beae4b4582e440555a4a603fb088a48a1906f4`.
- `scripts/adaptive_minimal_recovery_compare.py` — SHA256 `bd46578749b9d9b97f17cc730bb3adc60fc52d0af6a79dc958b97e435ca673d1`.
- `workspaces/adaptive_minimal_pool_20260923/REFERENCE_v1.json` — SHA256 `1f6dd2875c5c87ea3355f18c23c9e5b6f055e3cbb5ea12223727e16f2f523aba`.
- `workspaces/adaptive_minimal_pool_20260923/COMPARISON_v1.json` — SHA256 `65471e56b2897fed9e0d0206e7695247d3ba03e98b743c48ed59a005d22c95f4`.
- `workspaces/adaptive_minimal_pool_20260923/COMPARISON_recovered_v1.json` — SHA256 `47c30d3f4b50a0fa9da1db751900d3604d57b1648201214c404be7a2b361dafd`.
- `workspaces/adaptive_origin_recovery_20260923/run_v2/manifest.json` — SHA256 `d3f162d1d5f8834ee3c3dc8c5c785976df6451ee94b6a9140819fe2af0da7b2e`.
- `workspaces/adaptive_origin_recovery_20260923/run_v2/pool/manifest.json` — SHA256 `416ca2b119797eb6950640849fadf715191da25db187c9da9070a118f8e4a504`.
- `workspaces/adaptive_origin_recovery_20260923/run_v2/pool/final_collection.json` — SHA256 `e10ec1e158fdc30311f7ca65a01b48abaa4bd67834e5e8c776334fdb6c6c431c`.
