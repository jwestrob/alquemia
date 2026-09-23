# Standalone static transfer: reliable execution, no classification gain

**The fixed-origin standalone backend completes every prepared source, but does
not improve single-fold classification over the native released or recovered
adaptive methods.** Keep the existing scorers. Preserve standalone xTB as a
numerically qualified derivative backend for separately designed work.

All 832 logical solvent cells completed: 828 fresh calls and four exact pilot
reuses. No new MACE, DFT or geometry optimization ran. All 225 declared folds
remain in the report; 17 are unavailable because their original preparation
failed. No scientific failure, score substitution, label change or band refit
occurred in this branch.

## Frozen-reference transfer results

Each method uses its own previously frozen reference. The standalone reference
was fitted only to the designated 25 canonical examples before this transfer;
it retained 25/25 canonical and 3/3 consumed crystal calls. These 225 folds are
structural repeats of those 25 known reference proteins, not independent
biological observations or a fresh blind panel.

| Method | Correct | Wrong | Inconclusive | Unavailable | Total |
|---|---:|---:|---:|---:|---:|
| Native released static | 203 | 2 | 2 | 18 | 225 |
| Native recovered adaptive | 204 | 1 | 1 | 19 | 225 |
| Archived DFT | 199 | 3 | 6 | 17 | 225 |
| Standalone static | **201** | **2** | **5** | **17** | **225** |

Coverage matters. On the 207 sources shared with native released scoring,
standalone gives 200 correct / 2 wrong / 5 inconclusive versus 203 / 2 / 2.
On the 206 shared with recovered native adaptive, it gives 199 / 2 / 5 versus
204 / 1 / 1. Against DFT on all 208 prepared sources, standalone gives 201 / 2 / 5
versus 199 / 3 / 6. DFT uses its own core, Hamiltonian and reference; it is a
separate comparator, not the quantum component of the new score.

The sole older native solvent failure, A8R3S4 Ca-conditioned sample 3, is now
available and correctly Ca-supported (R = −405471.29857969936 model kcal/mol).
Its original failure remains in the historical row. This is useful coverage
recovery, not evidence of improved decisions on shared coverage.

## Actual noncorrect standalone calls

All names below refer to the existing source IDs, with seed −1. The raw values
and every old-method result remain in the full comparison.

| Source | Expected | Standalone outcome | R, model kcal/mol |
|---|---|---|---:|
| A0A3F2YLY8, Ca sample 1 | La | inconclusive | −405458.71831151 |
| A0A3F2YLY8, Ca sample 3 | La | inconclusive | −405457.68107656 |
| A0ACD6B9F2, Ca sample 4 | La | wrong Ca | −405475.65317953 |
| A0ACD6B9F2, La sample 4 | La | wrong Ca | −405466.34675777 |
| MMOL1770, Ca sample 1 | La | inconclusive | −405460.79914397 |
| A8R3S4, Ca sample 1 | Ca | inconclusive | −405465.99373974 |
| Q9Z4J7, La sample 4 | Ca | inconclusive | −405463.12212366 |

A0A3F2YLY8 Ca sample 1 improves from native-released wrong to inconclusive.
A0ACD6B9F2 La sample 4 changes from native-released inconclusive to wrong.
MMOL1770 Ca sample 1, A8R3S4 Ca sample 1 and Q9Z4J7 La sample 4 change from
native-released correct to inconclusive. No threshold adjustment rescues them.
Both wrong calls are structural samples of one protein family.

## Predeclared grouped summaries

Strict means every required source must be available. The summaries use the
existing arithmetic averaging rule and the same frozen method-specific bands.
They do not select the best fold or omit failed members.

| Standalone summary | Correct | Wrong / inconclusive | Unavailable | Total |
|---|---:|---:|---:|---:|
| La-conditioned four-fold mean | 23 | 0 / 0 | 2 | 25 |
| Ca-conditioned five-fold mean | 22 | 0 / 0 | 3 | 25 |
| Equal mean of the two arms | 21 | 0 / 0 | 4 | 25 |
| Every La three-of-four mean | 94 | 0 / 0 | 6 | 100 |

All available grouped decisions are correct. The 100 overlapping triples are
robustness checks, not 100 independent proteins. On shared grouped coverage,
standalone matches native released and adaptive decisions; extra availability
comes from recovering the A8R3S4 source. It changes the DFT Ca5 inconclusive
group to correct, while the other available DFT groups already agree.

## Method and numerical scope

For each fixed, exactly paired source geometry:

`E_M = native OMOL vacuum E_M + standalone GFN2 ALPB E_M − standalone GFN2 vacuum E_M`

`R = E_Ca − E_La`

Use standalone xTB 6.7.1, accuracy 0.02, fresh/no restart, electronic 300 K,
maximum 500 SCC iterations, exact charges and zero unpaired electrons. Both
solvent terms come from this executable. Pinned ALPB water settings retain
solvent 298.15 K, gsolv reference, P16/GBOBC, normal 230-point surface, hydrogen
bond/refshift terms and no ionic screening. This is not interchangeable with
the native ORCA GFN2 implementation by name alone. Raw model contrasts are
not binding free energies and have no universal zero.

Static frozen bands remain Ca maximum −405466.16155897034 and La minimum
−405456.4020100691 model kcal/mol. The failed standalone adaptive calibration
remains separate and closed; no optimized standalone pool is rescored here.
The original geometry, protonation, waters, source mappings and native MACE
receipts were reused exactly. All runtime parameter/state checks pass, with
24–108 SCC iterations and maximum charge closure error 1.30e−7 e.

The earlier 224-call experiment demonstrated energy-accuracy stability and
analytic/finite-difference force agreement on its declared development cases.
Those checks support this backend; they do not prove physical correctness or
force accuracy on every fold. This transfer adds coverage and discrimination
evidence without repeating a numerical sweep.

## Cost, tests and disposition

Job **1210508** completed in **244 s on 64 CPUs: 15,616 allocated core-seconds**,
zero GPU. Requested host memory was 128 GiB; scheduler batch MaxRSS was
781,304 KiB. Eight eight-thread molecular workers took 207.845 s collectively
(executor allocation 13,302.108 core-seconds); per-call wall median 1.800 s,
range 0.742–5.523 s. Startup, validation and collection account for the remaining
allocation time. Local preparation/report work and historical reused MACE/solvent
costs are additional. This is batch throughput, not a measured cold scanner time.

Four real-artifact tests pass in 10.500 s, zero skips. They cover exact
225/208/832/4/828 membership and reuse, preserved missing-source status, frozen
reference integrity (including a corrupted real-record rejection), actual
comparison counts and strict missing-member denominators. The result test was
explicitly unrun before molecular execution. No successful data were fabricated.

**Recommendation:** close standalone static as a replacement candidate under
this comparison. Keep its fast, working gradient backend for a separately
justified solvent-aware proposal experiment; retain native released and recovered
adaptive results. No production/default changes or further scoring are made.

Actual outputs: `workspaces/standalone_xtb_transfer_20260923/run_v1/` contains
`manifest.json`, `collection_v1.json`, `COMPARISON_v1.json`, `EXECUTION.json` and
`COSTS.json`. [Commands](COMMANDS.md), [frozen plan](PLAN.md) and [compact result](RESULT.json).
