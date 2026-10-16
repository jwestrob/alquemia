# LBT3 whole-scaffold saved comparison

Status: complete. Eight declared LBT3 cells; six saved Mex cells. No new model calls.

Raw energies are eV; differences below are kcal/mol. D = E(Dy) − E(La). Raw D is not an affinity.

| Source | Representation | La energy (eV) | Dy energy (eV) | D (kcal/mol) |
|---|---|---:|---:|---:|
| 7CCO | whole | -184236.161053 | -207459.369838 | -535539.916983 |
| 7CCO | matched_core | -36687.180673 | -59910.388048 | -535539.884447 |
| 7CCN | whole | -184234.051423 | -207460.384560 | -535611.966241 |
| 7CCN | matched_core | -36683.813882 | -59910.216753 | -535613.574328 |

| Balanced comparison | Δ (kcal/mol) |
|---|---:|
| 7CCO: D whole − matched core | -0.032535 |
| 7CCN: D whole − matched core | 1.608087 |
| whole: D 7CCN − 7CCO | -72.049259 |
| matched_core: D 7CCN − 7CCO | -73.689881 |

Rᵢ = D(Mex EFᵢ) − D(LBT3). Positive is conditional relative Mex La direction; no EF labels or summation.

| LBT3 source | Representation | R₁ | R₂ | R₃ |
|---|---|---:|---:|---:|
| 7CCO | whole | -3.776984 | 4.314385 | -3.801997 |
| 7CCO | matched_core | -3.809519 | 4.281850 | -3.834533 |
| 7CCN | whole | 68.272275 | 76.363644 | 68.247261 |
| 7CCN | matched_core | 69.880362 | 77.971731 | 69.855348 |

| Mathematical pool | La-selected source | Dy-selected source | D (kcal/mol) |
|---|---|---|---:|
| whole | 7CCO | 7CCN | -535563.317028 |
| matched_core | 7CCO | 7CCO | -535539.884447 |

Pools require matching chemical inventory and every cell; they are minima over two structures, not populations.
Whole/core changes include scaffold, D13/terminal charge and caps. Dry vacuum, terminal-state and physical-spin limitations remain. This consumed peptide control does not establish LanM profiling accuracy.

LBT3: 8 model calls, 22.480s application, job 12471; six historical Mex cells reused. Preparation's missing-H minimizations are recorded separately; their internal force-call count is unknown.
Exact inputs, runtime receipts, availability reasons and all matrices: `/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_coupled_response_20261016/lbt3_comparison_v2/execution/RESULT.json`.
