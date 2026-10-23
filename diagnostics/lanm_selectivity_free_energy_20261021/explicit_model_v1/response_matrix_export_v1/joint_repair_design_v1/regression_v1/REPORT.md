# Frozen water candidate: consumed trimer regression

**Decision: hold fresh validation and investigate collective interaction errors.**

The six-call diagnostic completes all eight metal/two-water cycles. Independent raw-token arithmetic reproduces all cycles and four Dy-minus-La contrasts. The frozen fit passes all25 development energy/force checks, yet7/8 collective energy errors grow. All8 individual dipole-error norms improve, while all4 paired dipole-error norms worsen.

| Cluster / water pair | Earlier paired energy error | Frozen candidate error |
| --- | ---: | ---: |
| W8 near | -1.4071 | -5.3119 |
| W8 far | -1.0393 | -4.6977 |
| W9 near | -2.2043 | -5.7713 |
| W9 far | -0.8789 | -4.5501 |

Errors are model minus the existing electronic reference, in kcal/mol. Every candidate paired three-body contrast is negative; every corresponding reference contrast is positive. This is evidence against treating the improved fitted totals as a repaired selectivity mechanism. It does not establish which functional form or parameterization will succeed.

No new physical threshold was introduced for these descriptive cycles. Preserve the candidate and all failures; no alternate iterate is selected. Fresh radial geometries and outcomes remain unconsumed. This diagnostic was released using the completed selected46 raw review and verified selection history. The full historical fit-receipt audit subsequently passed (`8658d5af`, root `1fe38e62`), with both strict-review failures retained. The scientific hold is unchanged.

Cost: six new analyze calls, Slurm18913 completed8seconds on oneCPU (8 allocated CPU-seconds). Independent saved arithmetic cost0.419CPU-seconds; no further molecular calls. Reference method/electronic-state limitations remain; these are finite consumed configurations, not a protein or aqueous exchange free energy.
