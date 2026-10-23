# Native component tradeoff: gain confirmed, one physical gate fails

The original Hamiltonian realizes nearly all of the predicted component gain when the four limits preserving earlier fitted scores are removed. This is evidence that the calibration policy restricted improvement. It does not establish a qualified model or require a new physical term.

| Quantity | Original local origin | Actual tested point |
|---|---:|---:|
| Pair component RMS error, kcal/mol | 26.6824 | 17.3193 |
| Collective remainder RMS error, kcal/mol | 22.3032 | 14.8343 |
| Pair Cartesian dipole RMS error, D | 0.614949 | 0.619410 |
| Collective Cartesian dipole RMS error, D | 0.0723095 | 0.0733691 |
| Susceptibility RMS error, a0^3 | 0.0341615 | 0.2083914 |

Actual component score q is0.4423864, close to the predicted0.4387129; the largest component prediction discrepancy is0.075454kcal/mol. All15force checks and9of10energy checks pass. The water-minus balanced-exchange error is−1.002200211kcal/mol and fails the unchanged±1limit. The small failure remains a failure; this point is not selected for sampling. Its inherited local prediction lay essentially on the boundary, so exact nonlinear acceptance is needed in any new calibration.

All four earlier block scores worsen (EF3.0303x, pair dipole1.0146x, collective dipole1.0295x, susceptibility37.2122x). Those origin-score ceilings were development safeguards, not physical tolerances. No new electrical accuracy threshold is imposed. Absolute electrical errors and all reference limitations remain visible.

This finite pilot is closed after exactly one vector,46cells and56native programs, including10analytic gradients and4response exports. Independent raw replay verifies56receipts/8matrix snapshots and maximum12overlapping programs. It cost11.913wall seconds and37.3597total CPU-seconds on the authorized login resource; independent saved review adds5.865wall/5.310CPU-seconds. No retry, optimizer, new reference or heldout case was used.

The component targets are paired La/Dy contrasts on consumed fixed water clusters, with provisional PBE0/ECP/common-basis and electronic-state limitations. The collective remainder contains all nonpair orders; it is not a unique polarization or charge-transfer measurement. Errors remain substantial and worse than the original pre-water-fit component baseline (~13.17/13.09kcal/mol). These are not exchange free energies.

Next: define one component-aware calibration with unchanged physical thresholds and nonlinear feasibility checks; report the electrical tradeoffs explicitly. The completed eight-step preservation-constrained fit remains closed. Periodic operator/stability development proceeds independently, with no dynamics admission from this pilot.

Evidence: native RESULTc361c162; independent reviewd0c0feb0; prospective policy9e994951; root actual review adjacent in the workspace. Original failed and completed evidence is unchanged.
