# Broad-coordinate response under unchanged per-system energies

6/6 finite searches admitted; 0/6 satisfy full-force convergence.

| System/source | Metal | Additional work (kcal/mol) | Final max/RMS force | Outcome |
| --- | --- | ---: | --- | --- |
| LBT 7CCO | La | -98.572405 | 2.141/0.921 | accepted_step_limit |
| LBT 7CCO | Dy | -104.193162 | 5.176/1.143 | accepted_step_limit |
| LBT 7CCN | La | -137.647592 | 9.740/1.988 | accepted_step_limit |
| LBT 7CCN | Dy | -112.937903 | 10.175/1.800 | accepted_step_limit |
| Mex 8FNS | La | -2945.113073 | 11.054/1.418 | accepted_step_limit |
| Mex 8FNS | Dy | -2942.034643 | 21.892/1.772 | accepted_step_limit |

Forces are kcal/mol/Å. Finite energy lowering is not a stationary minimum or accuracy test.

## Conditional exchange retained for continuity

**The scaffold Hamiltonians remain unmatched.** LBT uses whole-peptide MACE; Mex uses local electronics plus MM. These are not qualified preference predictions.

| LBT source | Old R | New R | Change (kcal/mol) |
| --- | ---: | ---: | ---: |
| 7CCO | 0.982619 | 9.681807 | 8.699188 |
| 7CCN | 25.103610 | 3.472352 | -21.631258 |

7CCN is Lu-conditioned. Mex retains three Nd spectators. No cross-source minimum, affinity estimate or experimental validation is implied.

## Decision

Use the measured component work and residual forces to diagnose missing mechanics. Module 3 must test matched local-electronic/MM representation; independent references must test accuracy. Energy lowering, a smaller source gap or a desired sign alone is insufficient.

Known scheduler allocation cost: 20904 CPU-s; missing receipts: none. Component calls, measured wall time and missing counters are preserved separately in COSTS.json.
