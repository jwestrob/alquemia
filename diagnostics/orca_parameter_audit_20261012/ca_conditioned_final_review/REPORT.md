# Corrected Ca-conditioned transfer — final review

All 125 declared sources are retained. Frozen canonical bands are unchanged; no new calibration or triples were created.

| Population / policy | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Full125 / origin | 104 | 1 | 4 | 16 |
| Full125 / operational | 105 | 1 | 1 | 18 |
| Full125 / mathematical | 105 | 1 | 1 | 18 |
| Common107 / origin | 102 | 1 | 4 | 0 |
| Common107 / operational | 105 | 1 | 1 | 0 |
| Common107 / mathematical | 105 | 1 | 1 | 0 |

Matched transitions: {"correct -> correct": 101, "correct -> inconclusive": 1, "inconclusive -> correct": 4, "wrong -> wrong": 1}.

## Remaining errors and abstentions

- `a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4`: wrong → wrong.
- `mmol_1770-pqq-la_model__conditioned_Ca__seed-1_sample-0`: correct → inconclusive.

All exclusions, interrupted native searches and changed decisions are listed in EXCEPTIONS.tsv; ALL125.tsv retains the full denominator. A missing accommodated pool is never filled with its valid origin score.

## Prior evidence remains separate

The existing 25-member canonical panel supplied these bands. Its fidelity is calibration evidence, not new predictive accuracy. Existing La-conditioned development rows and the established100 correlated triples are carried separately in RESULT.json; no new aggregate rule is invented.

## Components and measured cost

The audit reparsed actual backend outputs and checked unchanged component-v2 solvent differences, metal states, row minima, origin retention tolerance and all decisions. Historical absolute-v1 parity failures remain failed; this review does not reinterpret them.

{"absolute_v1_gate_failures_retained": 1300, "absolute_v1_maximum_error_Eh": 0.00019067857670052035, "components": 650, "maximum_component_error_Eh": 1.3395725773079903e-07, "maximum_population_error_e": 1.0496988167574273e-07, "unique_inputs": 1300}

Fresh solvent: 1264 xTB calls and 1264 tblite checks; 88725.616 measured process CPU-seconds. Stage wall time through collection: 7603.181 s; 91238.178 allocated CPU-seconds. These are different measures.

MACE recorded 3461 forwards and 263376.255 inference CPU-seconds. Interrupted forwards without terminal receipts retain unknown cost. No CPU-to-GPU extrapolation, end-to-end speedup, production promotion, or independent biological validation is claimed.

Recommendation: retain this corrected candidate as development evidence while keeping the remaining wrong/inconclusive cases and source dependence visible. Recalibration alone does not establish broader metal specificity.
