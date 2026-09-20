# Charged-context quantum check complete — 2026-09-20

All8native r2SCAN-3c/CPCM endpoints completed normally (job1202478).
The three charged Ca-family expansions shift R by only+1.11 to+3.28kcal/mol,
versus+49.13 to+51.03model-kcal in native vacuum OMOL. The charged La-family
expansion shifts+1.92 versus+19.78. No labels, coordinates or bands changed.

| Case | Class | Added charge | Native DFT context−core R | Native MACE context−core R | Included CPCM diagnostic change |
|---|---|---:|---:|---:|---:|
| a0acd6b9f2-pqq-la_model | La | -1 | 1.916343 | 19.776419 | -22.021967 |
| a8r3s4-pqq-la_model | Ca | -1 | 2.242472 | 51.025181 | -47.230373 |
| q88jh5-pqq-la_model | Ca | -1 | 1.111000 | 49.128656 | -39.141605 |
| q9z4j7-pqq-la_model | Ca | -1 | 3.278275 | 49.391381 | -35.897893 |
| 1H4I | Ca | 0 | -2.032777 | -22.850199 | 32.476439 |
| 4MAE | La | 0 | -4.425467 | -35.247174 | 30.834959 |

Converted energy contrasts use kcal/mol for DFT and model-kcal for MACE.
The CPCM column is already part of the native DFT total, never added again.
Different densities and methods prevent interpreting its subtraction as an actual
vacuum-DFT calculation or assigning the discrepancy uniquely to solvent. The
opposing solvent response nonetheless motivates the separately approved matched
GFN2ALPB−vacuum transfer for vacuum MACE.

Cost:3415s×64allocatedCPUs=218560core-s;0GPU-s. Local observer cost1.823CPU-s
is separate. All outputs, execution receipts and collection are pinned in
COMPLETION.json and collected_1202478.json. Root read these actual results before
proposing the next experiment; no repeated DFT or numerical rescue was launched.

Next approved scope: ../compact_solvation_20260920/PLAN.md.
