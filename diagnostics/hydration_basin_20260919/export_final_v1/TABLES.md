# Coupled water response: actual native checks

All energy changes are kcal/mol. A negative change lowers the native DFT energy. These are electronic response components; occupancy remains unavailable.

| Round | Center | Native change this step | Total change | Prediction error | Energy check | Native gradient max |
|---:|---|---:|---:|---:|---|---:|
| 0 | 1F6S__11__Ca | -4.20764 | -4.20764 | -0.01901 | pass | 3.2356 |
| 0 | 1F6S__11__La | -2.46261 | -2.46261 | 0.04894 | pass | 2.1542 |
| 0 | 6IP9__110__Ca | -2.05571 | -2.05571 | 0.07029 | pass | 1.2815 |
| 0 | 6IP9__110__La | -1.81313 | -1.81313 | 0.05784 | pass | 0.6114 |
| 1 | 1F6S__11__Ca | -0.34096 | -4.54859 | 0.02296 | pass | 0.7264 |
| 1 | 1F6S__11__La | -0.15204 | -2.61464 | 0.02057 | pass | 0.6344 |
| 1 | 6IP9__110__Ca | -0.25486 | -2.31057 | 0.02180 | pass | 0.3890 |
| 1 | 6IP9__110__La | -0.05297 | -1.86610 | 0.01403 | pass | 0.2624 |
| 2 | 1F6S__11__Ca | -0.06069 | -4.60929 | 0.00840 | pass | 0.7024 |
| 2 | 1F6S__11__La | -0.01780 | -2.63244 | 0.00098 | pass | 0.2815 |
| 2 | 6IP9__110__Ca | -0.01616 | -2.32673 | 0.00532 | pass | 0.1737 |

Gradient maxima use the physical coordinates: kcal/mol/Å for translations, kcal/mol/radian for rotations; the reported maximum uses the declared 1 Å/radian numerical scale.

| Round | Site | Total response contribution to Ca−La | Paired step prediction error |
|---:|---|---:|---:|
| 0 | 1F6S 11 | -1.74503 | -0.06795 |
| 0 | 6IP9 110 | -0.24258 | 0.01245 |
| 1 | 1F6S 11 | -1.93395 | 0.00239 |
| 1 | 6IP9 110 | -0.44447 | 0.00777 |
| 2 | 1F6S 11 | -1.97685 | 0.00743 |
| 2 | 6IP9 110 | unavailable | unavailable |

Larger Ca−La is more La-like within this protocol. The aquo offset cancels from these response contributions; historical decision bands do not transfer.

## Latest measured geometries, with explicit endpoint reuse

| Site | Ca round | La round | Total Ca−La response contribution | Both checked interior minima |
|---|---:|---:|---:|---|
| 1F6S 11 | 2 | 2 | -1.97685 | no |
| 6IP9 110 | 2 | 1 | -0.46063 | yes |

These use actual native energies. Earlier qualified endpoints are reused explicitly when only their partner receives another step. There is no baseline substitution.

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
| 1 | 1F6S__11__Ca | 0.2387 | 0.3147 | no |
| 1 | 1F6S__11__La | 0.2395 | 0.2928 | no |
| 1 | 6IP9__110__Ca | 0.3426 | 0.4068 | no |
| 1 | 6IP9__110__La | 0.2796 | 0.3699 | no |
| 2 | 1F6S__11__Ca | 0.2781 | 0.3609 | no |
| 2 | 1F6S__11__La | 0.2413 | 0.2964 | no |
| 2 | 6IP9__110__Ca | 0.3502 | 0.3988 | no |

The domain remains 0.20 Å / 0.35 rad. These are conditional harmonic-model diagnostics, not validated physical fluctuations or entropy estimates.

## Occupancy remains unavailable

- validated stationary coupled basin spectra and their physical extent.
- bound internal water vibration and its consistent reference.
- non-electrostatic solvent contribution.
- competing orientation basins and all occupancy-state corrections.

No new biological accuracy claim, default change, threshold fit or PQQ rescore.
