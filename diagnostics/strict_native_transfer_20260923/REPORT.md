# Full225 strict-native scalar transfer

**207 correct, 0 wrong, 1 inconclusive and 17 unavailable.** The method uses the same completed tenfold-precision geometry pools and MACE energies, with fresh strict native solvent calculations. Production is unchanged.

| Method | Correct / wrong / inconclusive / unavailable |
|---|---|
| Released static | 203 / 2 / 2 / 18 |
| Prior precision tenfold | 206 / 0 / 2 / 17 |
| Strict fresh tenfold | 207 / 0 / 1 / 17 |

## Matched coverage and conditioning

On the same 207 available sources, Released static: 203 / 2 / 2 / 0; strict: 206 / 0 / 1 / 0.
On the same 208 available sources, Prior precision tenfold: 206 / 0 / 2 / 0; strict: 207 / 0 / 1 / 0.
La100: strict 97 / 0 / 1 / 2.
Ca125: strict 110 / 0 / 0 / 15.

| Strict aggregate | Released | Prior precision | Strict fresh |
|---|---|---|---|
| La4 | 23 / 0 / 0 / 2 | 23 / 0 / 0 / 2 | 23 / 0 / 0 / 2 |
| Ca5 | 21 / 0 / 0 / 4 | 22 / 0 / 0 / 3 | 22 / 0 / 0 / 3 |
| balanced | 20 / 0 / 0 / 5 | 21 / 0 / 0 / 4 | 21 / 0 / 0 / 4 |
| La100_triples | 94 / 0 / 0 / 6 | 94 / 0 / 0 / 6 | 94 / 0 / 0 / 6 |

Each aggregate requires every declared member. The100 triples and225 folds are correlated samples of25 consumed proteins, not independent biological validations. All17 original preparation exclusions remain; no failed cell is filled with an old energy.

## Decision changes versus prior precision

- `q4w6g0-pqq-la_model__conditioned_Ca__seed-1_sample-2`: inconclusive → Ca-supported; ΔR=-7.987209876 model kcal/mol; new score under old bands: Ca-supported.

Remaining wrong/inconclusive calls:

- `c5axv8-pqq-la_model__conditioned_La__seed-1_sample-3`: inconclusive (La-class reference), R=-405457.848334666.

## Numerical interpretation

All geometry and MACE cells were reused unchanged. 2469/2496 scalar cells remain within0.1 kcal/mol of their original loose result; 201/208 pools remain within0.2. Maximum absolute component/pool shifts: 44.512340636/7.987209876 kcal/mol. These comparisons describe the change; reproducing the loose solution is not an acceptance criterion.

Among45 complete conditioning groups, spread shrinks in30, grows in15 and is unchanged in0. Median range: 10.822193→10.760895 model kcal/mol. Whole-group medians, raw contrasts, selected candidates and each component remain in the full JSON.

The frozen fresh strict32 reference uses only the original25 calibration sources. No threshold was refitted here. Strict32 fresh/seeded agreement is a separate numerical qualification; these transfer sources receive one fresh start each. Scalar consistency does not qualify forces, prove a unique electronic solution or turn electronic contrasts into binding free energies.

## Actual cost and reproducibility

Four CPU-only32-core/64GiB jobs: 1211337, 1211338, 1211339, 1211340. New calls: 2400/2400 complete; eight whole pools/96 cells reused separately. No new MACE, searches, DFT or preparation chemistry. Allocated 168352 core-seconds; actual accounting steps 18960.007 CPU-seconds; zero requested GPU-seconds. Summed native ORCA time 147909.549s. Full job/step receipts, failures and per-cell times are in COSTS.json.

Local source/manifest audit took150.387s; implementation, tests and reporting add unmetered local work. Reused qualification/comparator allocations and historical MACE work are excluded from these new-job totals, not declared free. This is not a matched hardware speed comparison.

Full source/aggregate results: `workspaces/strict_native_transfer_20260923/run_v1/COMPARISON.json`; compact summary, costs and pinned receipts beside it. See COMMANDS.md for the exact no-molecule comparison operation.
