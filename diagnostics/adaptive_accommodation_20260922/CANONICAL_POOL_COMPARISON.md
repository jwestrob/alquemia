# Adaptive canonical-pool comparison

**The four-mode adaptive proposals produce real energy changes, but do not yet demonstrate better discrimination.** The full 30-source comparison has 29 available pairs. One canonical member fails during optimization, so the strict adaptive reference remains unavailable; no 24-member calibration was fitted.

## Outcome under unchanged static bands

These are **transfer checks**, not calibrated adaptive affinity or class predictions. The same original 25 canonical members remain in the denominator.

| Population | Static context, same available sources | Earlier shared pool, same available sources | Adaptive, same available sources | Adaptive unavailable |
|---|---|---|---|---|
| Canonical25 (24 common) | 24 correct | 23 correct / 1 inconclusive | 22 correct / 1 wrong / 1 inconclusive | 1 |
| Crystal3 (3 common) | 3 correct | 3 correct | 3 correct | 0 |
| Unknown PLM2 | labels unknown | Ca / inconclusive | Ca / inconclusive | 0 |

- **Q9Z4J7** (Ca reference): earlier shared-pool inconclusive → adaptive La-supported under static bands; R changes −405462.023718→−405459.009899.
- **Q88JH5** (Ca reference): earlier shared-pool Ca-supported → adaptive inconclusive; R changes −405466.838456→−405459.555871.
- **MMOL1770**: La optimizer attempted an unsupported trial geometry and returned no candidate. The valid original coordinates replayed exactly; this is an optimization failure, not evidence of an invalid source structure. Its case remains unavailable, including the otherwise successful Ca proposal.
- The unknown PLM contrasts move +2.516051 and +0.305201 model kcal/mol from the prior shared pool. Their old-band calls remain Ca-supported and inconclusive; neither supplies a correctness label.
- Mathematical and operational variants give identical old-band calls. Their raw values remain separate (1KB0 differs slightly), and both new-reference decisions remain unavailable for every case.

## Reference and numerical accounting

- Protocol: `nikasha_adaptive_angular_common_geometry_native_OMOL_GFN2_ALPB_v1`.
- Separate strict reference artifact: `remaining26_pool_v1/REFERENCE.json`, reference ID `Nikasha_adaptive_angular_common_geometry_canonical25_reference_v1`. Both variants have status `unavailable_calibration_member`, 24/25 available, and null bands, extrema and gap. No missing member was removed or baseline score substituted.
- The actual pilot retains its original 31/32 converged solvent targets and failed 4MAE vacuum output. The explicit MaxIter500 diagnostic reproduced two controls and completed the missing target without changing parameters/convergence tolerances; its three attempts and primary failed case remain pinned.
- The remaining 26 extension uses the documented same-Hamiltonian MaxIter500 policy and the same qualification pin. All 200 new GFN2 targets and 50 cross-MACE calls complete. Its proposal run returned 51/52 endpoint candidates, yielding 25/26 complete paired pools.
- Across pilot and continuation: 60 attempted endpoint optimizations, 59 candidates; 58 cross-MACE evaluations; 235 actual GFN2 attempts (234 successful, one retained primary failure), representing 232 completed target cells plus two control repetitions. No DFT, new fold, chemical-state change or default change.

## Measured cost and reusable outputs

- Full adaptive continuation including pilot proposals, failed solvent attempt and numerical diagnostic: **60,888 allocated core-seconds / 376 GPU-seconds**.
- Remaining26 matrix completion alone: **42,016 core-seconds / 27 GPU-seconds**. These costs exclude original preparation, origins, earlier restricted proposals and base-pool evaluation; mixed CPU hosts do not support a matched-hardware speed claim.
- Final collection: `workspaces/adaptive_accommodation_20260922/remaining26_pool_v1/final_collection.json`.
- Every raw pair/old-band transfer: `original30_comparison.json` and `.md`; direct static/base/adaptive columns: `paired_contrasts.csv`; counts/costs: `summary.json` and `scheduler_accounting_final.txt` in the same directory.
- Both primary/recovery metadata and original labels remain traceable. Thirteen real-artifact comparison tests pass, including the actual 30-source null-reference path; no fabricated scientific fixtures.

**Interpretation:** four angular modes change computed accommodation energies, but these consumed reference results do not establish discriminatory improvement. The common pool is technically working; the optimizer has one explicit coverage failure, and the old-band class transfer worsens for two Ca controls. Keep the released model unchanged. Parent will decide any separately declared next mechanical test; no follow-on was launched here.
