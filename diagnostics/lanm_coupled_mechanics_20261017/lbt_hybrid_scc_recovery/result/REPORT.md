# Isolated-core numerical recovery passes

All four declared GFN2 calls and three energy/gradient comparisons pass. The original failed pool remains unavailable until the replacement component is explicitly assembled and all37 untouched cells are evaluated.

| Comparison | Energy difference Eh | Max physical gradient difference kcal/mol/A | Pass |
| --- | ---: | ---: | --- |
| target_branch_consistency | 9.66338121e-13 | 2.95032276e-06 | True |
| preceding_control | 9.66338121e-13 | 2.89696082e-06 | True |
| origin_control | 9.66338121e-13 | 3.55271985e-06 | True |

Fixed limits:1e−7Eh and0.005kcal/mol/A; original native SCC/state/parameter/gradient checks unchanged. All actual initial damping settings are printed by the native executable. Failedpoint0.05/0.10 converges in112/157iterations; controls263/59. Only the already declared0.05 target component is nominated for recovery. This is not a universal fallback or evidence of a unique physical electronic state.

Cost:4GFN calls,54.708061application seconds,52.440934native CPU-s; scheduler55s on1CPU. No MACE,MM,ALPB,DFT or optimization. Preflight13008 is separate.

Original12992 isolated SCC failure and all raw outputs remain immutable. A new declared component assembly must preserve its accepted MACE/embeddedGFN/MM and exact ALPB cache, with explicit solver-exception provenance. The remaining37cells keep the original evaluator and stop-on-failure rule; no repeat of66successful cells.
