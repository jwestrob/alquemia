# Original 30 shared-pool result

**Giving both metals the same three geometries changes small energy contrasts but adds no classification improvement on these consumed references.** All 30 sources remain available; no molecular calculation failed. This finite-pool fairness check preserves the prior broad separation, and the separate fold comparison is still required. Those folds repeat the same biological groups.

## Scientific result

- The shared set is original geometry plus the existing Ca- and La-proposed geometries. Both metals receive the same set; no new optimization or DFT was performed.
- With unchanged original static bands: 24/25 canonical calls correct, one inconclusive (Q9Z4J7), and 3/3 crystal calls correct. These are unchanged from own-proposal selection.
- Separate canonical-only recalibration gives 25/25 by construction and retains 3/3 crystal calls. This is developmental calibration, not fresh validation.
- Both mathematical and operational references have Ca_max **−405462.0237183686** and La_min **−405456.53376105055** model kcal/mol, gap **5.4899573181**. The prior own-proposal gap was 5.8443969721; the original static gap was 5.0763662803.
- Both unknown PLM predictions remain unchanged: 83440678cbbd658047c9 Ca-supported; 07ab500e3df76b30d71c inconclusive. They supply no truth labels.

Only the following raw contrasts differ from own-proposal selection; all other cases are identical:

| Case | Mathematical ΔR | Operational ΔR |
|---|---:|---:|
| 1KB0 | -0.071212972 | -0.071212972 |
| a0acd6b9f2-pqq-la_model | +0.093666266 | +0.000000000 |
| c5axv8-pqq-la_model | -0.354439654 | -0.354439654 |
| p12293-pqq-la_model | -0.087250969 | +0.000000000 |
| q4w6g0-pqq-la_model | +0.116061716 | +0.116061716 |
| q60ar6-pqq-la_model | -0.032855054 | +0.000000000 |
| q8gr64-pqq-la_model | +0.178925768 | +0.178925768 |

ΔR is model kcal/mol; positive is more La-like. Mathematical row minima include tiny improvements, while the operational rule requires >0.10 kcal/mol improvement over the original geometry. The two variants have separate frozen reference IDs even though their extrema coincide here.

## Executed work and cost

- Pilot 4: GPU 1209840 (8 cross-MACE) and CPU 1209845 (16 GFN2). Remaining 26: GPU 1209850 (52 cross-MACE) and CPU 1209862–65 (104 GFN2). All 120 GFN2 attempts converged and terminated normally; all 60 MACE calls succeeded.
- Pending 1209841/1209851–54 were canceled before execution; zero elapsed allocation and no molecular attempts. Existing submission receipts retain the placement changes.
- Incremental allocation: **25,472 core-seconds and 38 GPU-seconds**. GPU drivers measured 32.728557 s. This excludes earlier geometry searches, reused origins and preparation; it is not the end-to-end production cost.
- Final collection/calibration/inspection performed no molecular calls. Production and original bands remain unchanged.
- Seven actual-artifact comparison tests pass, including the completed-pool test; zero skips.

## Immutable outputs and next comparison

- Reference: `workspaces/nikasha_shared_pool_20260922/remaining26_v1/REFERENCE.json`; SHA256 `f96b2425202c0c24f11312f28109286cb88b5fb520a75a16d2a9e4ee9dbd0911`; frozen `2026-09-22T22:56:26.688183+00:00` before full-fold execution.
- Actual remaining 26: `workspaces/nikasha_shared_pool_20260922/remaining26_v1/final_collection.json` (26/26 available, 52/52 MACE, 104/104 GFN2).
- Every 30-case result, choice, old/new call: `workspaces/nikasha_shared_pool_20260922/remaining26_v1/original30_comparison.json` and `.md`; compact counts/costs: `original30_summary.json`.
- Compare the already-declared 225 folds with this frozen reference and `workspaces/nikasha_recovery_20260922/proposal_comparison.json`; do not recalibrate on folds. All sources, failures, conditioning strata, strict medians and 100 triples remain in scope.
