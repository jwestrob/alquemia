# Interpretation of the fixed block bounds

The declared six-matrix calculation is complete. Result SHA256: `a1493599d8cf9e15a3e04b6ecf2d104ee41e9590632cc01666b9b4c739c32e7d`, at the matching workspace `actual_v1/RESULT.json`. Source `2ce80ea7`, declaration `b20e2297`, inputs `9b176ba8` and READY `a607a9fc` were frozen before evaluation.

Both fixed partitions pass the unchanged sufficient criterion on all four one-/two-water fragments. Both fail on the two W9 centers, despite saved dense minimum eigenvalues of 0.611101 (La) and 0.611303 (Dy).

| W9 center | Site lower bound | Intact-water lower bound | Saved dense minimum |
|---|---:|---:|---:|
| La | −0.734264 | −0.230154 | 0.611101 |
| Dy | −0.587827 | −0.204869 | 0.611303 |

Water grouping improves the estimate but is already too conservative for these small vacuum clusters. In the water partition, the limiting block is the first intact water for both metals: its diagonal minimum is 0.619912, while the summed off-block norms are 0.850066 (La) and 0.824781 (Dy). The estimate replaces the directional structure of separate couplings with their summed worst-case magnitudes. Its failure therefore provides no evidence of an actual unstable mode.

These bounds may serve as an inexpensive sufficient screen where they pass. They cannot serve as the sole stability test for the admitted W9 systems, and this evaluation does not support extrapolating them to periodic solvent. A stronger certificate construction or retained fallback is needed; no new partition, weighting or threshold was tried here.

All six saved dense spectra replay within 3.34e−15. The largest symmetrization correction has operator norm 4.17e−17, far below the negative-bound deficits. These ordinary floating-point calculations are not rigorous rounded certificates. Cost: 0.559 wall seconds / 0.510 CPU seconds, one thread, zero native calls.
