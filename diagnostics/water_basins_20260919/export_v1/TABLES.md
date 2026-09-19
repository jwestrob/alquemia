# Finite water basins: actual sampling and native checks

These are conditional cheap-potential configurational integrals on declared common physical domains. They are not occupancy, absolute entropy, calibrated discrimination scores, or unbiased DFT reweighting.

| Center | Native representative weighted error | Maximum error | Native gate | MACE calls |
|---|---:|---:|---|---:|
| 1F6S__11__Ca | 0.56701 | 0.71975 | False | 1744 |
| 1F6S__11__La | 0.33096 | 0.73101 | False | 1882 |
| 6IP9__110__Ca | 0.42882 | 0.75950 | False | 1112 |
| 6IP9__110__La | 0.68470 | 2.15629 | False | 1516 |

| Center | Domain Å/rad | Conditional Fconfig | ESS | Scramble difference | Outer10% weight | Numerical gate |
|---|---|---:|---:|---:|---:|---|
| 1F6S__11__Ca | 0.3/0.5 | 8.93622 | 121.5 | 0.11582 | 0.2141 | False |
| 1F6S__11__Ca | 0.45/0.8 | 8.70477 | 220.2 | 0.04684 | 0.0522 | True |
| 1F6S__11__Ca | 0.6/1.1 | 8.61956 | 217.8 | 0.07947 | 0.0662 | True |
| 1F6S__11__La | 0.3/0.5 | 9.11680 | 134.1 | 0.05047 | 0.1387 | True |
| 1F6S__11__La | 0.45/0.8 | 8.87776 | 250.1 | 0.04205 | 0.0551 | True |
| 1F6S__11__La | 0.6/1.1 | 8.82466 | 257.7 | 0.04841 | 0.0177 | True |
| 6IP9__110__Ca | 0.3/0.5 | 7.41699 | 29.5 | 0.35617 | 0.3446 | False |
| 6IP9__110__Ca | 0.45/0.8 | 6.59027 | 85.3 | 0.00284 | 0.2369 | False |
| 6IP9__110__Ca | 0.6/1.1 | 6.26439 | 134.0 | 0.03765 | 0.1931 | True |
| 6IP9__110__La | 0.3/0.5 | 8.15152 | 62.0 | 0.09973 | 0.3034 | False |
| 6IP9__110__La | 0.45/0.8 | 7.58335 | 156.4 | 0.06064 | 0.1629 | True |
| 6IP9__110__La | 0.6/1.1 | 7.33903 | 115.4 | 0.05793 | 0.1482 | False |

| Pair | Domain Å/rad | Cheap Ca−La Fconfig contribution | Both numerical gates | Both native gates |
|---|---|---:|---|---|
| 1F6S__11 | 0.3/0.5 | -0.18059 | False | False |
| 1F6S__11 | 0.45/0.8 | -0.17298 | True | False |
| 1F6S__11 | 0.6/1.1 | -0.20510 | True | False |
| 6IP9__110 | 0.3/0.5 | -0.73453 | False | False |
| 6IP9__110 | 0.45/0.8 | -0.99307 | False | False |
| 6IP9__110 | 0.6/1.1 | -1.07464 | False | False |

All energies in kcal/mol. Configurational values use the explicit common measure in the agreement; its additive gauge cancels in these equal-water-count pairs. Kinetic/internal/quantum/nonpolar terms are not supplied. A failure of broad-domain target agreement does not negate the previously measured local relaxation improvement.
