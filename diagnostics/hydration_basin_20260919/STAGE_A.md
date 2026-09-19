# Coupled water calculation — Stage A complete

All 20 initial and proposed-minimum numerical curvature checks pass. Four cheap
minima are interior; 16 reach the translation boundary. No entropy or occupancy
is inferred from these results. Native validation is running.

| State | Negative initial modes | Translation (Å) | Rotation (rad) | Predicted ΔE (kcal/mol) | Interior minimum |
|---|---:|---:|---:|---:|:---:|
| 1F6S__01__Ca | 0 | 0.2000 | 0.0961 | -2.2124 | no |
| 1F6S__01__La | 0 | 0.1835 | 0.0349 | -0.6044 | yes |
| 1F6S__10__Ca | 0 | 0.2000 | 0.0510 | -2.3407 | no |
| 1F6S__10__La | 0 | 0.2000 | 0.1431 | -1.7225 | no |
| 1F6S__11__Ca | 0 | 0.2000 | 0.0926 | -4.2266 | no |
| 1F6S__11__La | 0 | 0.2000 | 0.1808 | -2.4137 | no |
| 6IP9__001__Ca | 0 | 0.2000 | 0.1491 | -1.3196 | no |
| 6IP9__001__La | 0 | 0.1316 | 0.1114 | -0.9921 | yes |
| 6IP9__010__Ca | 0 | 0.2000 | 0.1280 | -0.9591 | no |
| 6IP9__010__La | 0 | 0.2000 | 0.1338 | -0.4132 | no |
| 6IP9__011__Ca | 0 | 0.2000 | 0.2398 | -3.0782 | no |
| 6IP9__011__La | 0 | 0.2000 | 0.2051 | -4.0118 | no |
| 6IP9__100__Ca | 0 | 0.2000 | 0.0419 | -1.0096 | no |
| 6IP9__100__La | 0 | 0.1875 | 0.0456 | -1.1332 | yes |
| 6IP9__101__Ca | 0 | 0.2000 | 0.1477 | -2.4450 | no |
| 6IP9__101__La | 0 | 0.1824 | 0.1295 | -2.6772 | yes |
| 6IP9__110__Ca | 0 | 0.2000 | 0.1489 | -1.9854 | no |
| 6IP9__110__La | 0 | 0.2000 | 0.1435 | -1.7553 | no |
| 6IP9__111__Ca | 0 | 0.2000 | 0.2476 | -4.3933 | no |
| 6IP9__111__La | 0 | 0.2000 | 0.2034 | -6.2065 | no |

1966 objective evaluations; 408 GPU-seconds and 6528 allocated core-seconds.
These are model predictions, not native optimized energies or new biological tests.
