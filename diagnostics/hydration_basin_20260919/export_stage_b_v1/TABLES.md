# Coupled water response: actual native checks

All energy changes are kcal/mol. A negative change lowers the native DFT energy. These are electronic response components; occupancy remains unavailable.

| Round | Center | Native change this step | Total change | Prediction error | Energy check | Native gradient max |
|---:|---|---:|---:|---:|---|---:|
| 0 | 1F6S__11__Ca | -4.20764 | -4.20764 | -0.01901 | pass | 3.2356 |
| 0 | 1F6S__11__La | -2.46261 | -2.46261 | 0.04894 | pass | 2.1542 |
| 0 | 6IP9__110__Ca | -2.05571 | -2.05571 | 0.07029 | pass | 1.2815 |
| 0 | 6IP9__110__La | -1.81313 | -1.81313 | 0.05784 | pass | 0.6114 |

Gradient maxima use the physical coordinates: kcal/mol/Å for translations, kcal/mol/radian for rotations; the reported maximum uses the declared 1 Å/radian numerical scale.

| Round | Site | Total response contribution to Ca−La | Paired step prediction error |
|---:|---|---:|---:|
| 0 | 1F6S 11 | -1.74503 | -0.06795 |
| 0 | 6IP9 110 | -0.24258 | 0.01245 |

Larger Ca−La is more La-like within this protocol. The aquo offset cancels from these response contributions; historical decision bands do not transfer.

| Center | Direction | DFT even energy | Predicted even energy | Error | Frozen tolerance | Result |
|---|---|---:|---:|---:|---:|---|
| 1F6S__11__Ca | soft | 0.012884 | 0.003137 | -0.009748 | 0.005000 | fail |
| 1F6S__11__Ca | response | 0.108352 | 0.109159 | 0.000807 | 0.027088 | pass |
| 1F6S__11__La | soft | 0.007916 | 0.007882 | -0.000034 | 0.005000 | pass |
| 1F6S__11__La | response | 0.118064 | 0.127566 | 0.009502 | 0.029516 | pass |
| 6IP9__110__Ca | soft | 0.011469 | 0.008439 | -0.003029 | 0.005000 | pass |
| 6IP9__110__Ca | response | 0.096517 | 0.097360 | 0.000843 | 0.024129 | pass |
| 6IP9__110__La | soft | 0.012900 | 0.014257 | 0.001357 | 0.005000 | pass |
| 6IP9__110__La | response | 0.100340 | 0.099601 | -0.000738 | 0.025085 | pass |

## Conditional thermal extent of the cheap basin model

| Round | Center | Largest translation RMS (Å) | Largest rotation RMS (rad) | All RMS within local domain |
|---:|---|---:|---:|---|
| 0 | 1F6S__11__Ca | 0.2853 | 0.3411 | no |
| 0 | 1F6S__11__La | 0.2395 | 0.2832 | no |
| 0 | 6IP9__110__Ca | 0.3299 | 0.4338 | no |
| 0 | 6IP9__110__La | 0.2895 | 0.3948 | no |

The domain remains 0.20 Å / 0.35 rad. These are conditional harmonic-model diagnostics, not validated physical fluctuations or entropy estimates.

## Occupancy remains unavailable

- validated stationary coupled basin spectra and their physical extent.
- bound internal water vibration and its consistent reference.
- non-electrostatic solvent contribution.
- competing orientation basins and all occupancy-state corrections.

No new biological accuracy claim, default change, threshold fit or PQQ rescore.
