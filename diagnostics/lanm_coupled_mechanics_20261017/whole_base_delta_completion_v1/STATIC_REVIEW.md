# Independent review: one regional MACE forward and saved assembly

19 October 2026. **PASS; no static blocker found for the declared GPU stage and subsequent saved-only assembly.** This review inspected source, saved metadata and file hashes only. It performed no scientific imports, numerical execution or submissions.

Frozen artifacts:

- COMPLETE.py: `05414274fc4139d7e274f56f1bf513f76c20bed696397890a6eeee99b652dc2c`.
- SPEC.md: `6709cd34188189570bd3df5329df976cec4b3112d5809ebce80e6e1acfa068d4`.
- CONFIG.template.json: `07924cc0a6c4f3f2c400ef6e6950bbddd3bee16b6effd9690adf4adaf4aaeb5e`.
- Bound CONFIG_GPU.json: `114233afceeaa21b17d96498ae78b0f768cc8e619b9716935ff3b1b405d1b3b4`.

All exact path/SHA records in the bound config verify. Its only change from the template is the actual CPU result pin. Original preflight/scout BASE_DELTA and CPU sources remain unchanged.

## Prerequisite and call ceiling

Actual CPU job 14802 result `ba406175bcc4999757fd08210a1e52bbe583461b1cee5b92537a4042e6f0aceb` reports four accepted requests, status `development_cpu_complete`, cumulative five calls including the scout, and all four 7CCN slots not requested. The source checks that exact stage/source/config, counts and ordered request set before model loading. Saved CPU preparation replays the original eight states and scout; all five new GFN receipts are then re-admitted against native settings, coordinates, physical IDs and the regional map.

The GPU stage invokes only the frozen unbound `HybridEvaluator.mace` method, using the qualified NativeMACE constructor without initializing a hybrid evaluator. It requests one gradient forward for the missing baseline-172 7CCO after-Dy region: target 0, charge −2, physical multiplicity 6, exact archived geometry/caps, float64 and original model/head. It does not invoke the whole-256 evaluate method, a warmup forward, cached regional pairs, MM, embedded GFN or ALPB. Native forward and attempted-call counts must equal one. Cumulative development count is then six; four conditional whole-GFN challenge calls remain outside this executable, under the original ceiling of ten.

Failed preparation cannot load the model. Failed GPU execution retains the attempted count, raw receipt and explicit attempted/pre-forward unavailable row. All eight slots exist from initialization. No fallback, retry or automatic challenge is implemented. Root owns the one-GPU allocation and its separate launcher/cost record.

## Raw gradients and assembly contract

Collection invokes saved CPU preparation and native replay only; no model constructor or molecular engine is used. The three cached 7CCO local pairs are rebuilt from their raw MACE forces and isolated-GFN gradients by frozen BASE_DELTA.prepare. The missing pair uses the same independent conversion/scatter order: negative MACE force times EV, GFN gradient times HA divided by BOHR, both cap anchors, zero field gradient, and negation of the scattered regional GFN. No old full hybrid gradient is reused.

The full-256 order is GFN_full, MACE_region, negative_GFN_region, fixed_FF_ALPB. Candidate gradient is their left-to-right sum; base gradient is whole GFN plus ALPB; forces are negative gradients. Raw energies are retained and scalar totals use math.fsum. The saved comparator takes native before/after differences before conversion, rather than subtracting these large converted totals.

I compared the final emitted schema with frozen COMPARE.py `63f1f36a…`: root design versus adapter specification keys, protocol, stage `collect_7CCO`, status `complete_7CCO`, eight ordered rows, native energy keys, component receipts, state/map fields and exact component closure agree. Partial assembled files now retain their pins separately; the four complete rows publish only after all four assemble. The comparator therefore cannot produce survivor work from a failed parent. Missing before-Dy whole-MACE and after-Dy old-hybrid arrays remain missing.

## Remaining actual gates

This is source readiness, not GPU or assembly admission. The GPU result must exist and pass its original checks before binding collection; the completed assembly must then bind the independent comparator. Neither this review nor native admission satisfies the frozen both-metal 1 kcal/mol work triage, derivative qualification, physical-state validity or force accuracy. No optimizer or 7CCN challenge follows automatically.
