# Whole-peptide DFT force scout: admitted result

17 October 2026 UTC. Native job13051 and the original saved collector13283
completed. The complete256-atom before-La gradient passes the unchanged native
admission and exact saved replay. One reference cell is available; the other
three fixed before/after La/Dy cells have not run.

## Finding

The conditional PBE0 reference does not support the very large terminal force
produced by the hybrid model. Whole-peptide MACE is closer to this reference
over the full force array, although substantial differences remain.

All three composite models use the same fixed-charge ALPB term at the same
coordinates. Values below are kcal/mol/Angstrom, copied from the saved
`FORCE_COMPARISON.json`; projections use directions frozen before this DFT run.

| Quantity | Whole PBE0 + ALPB | Hybrid | Whole MACE + ALPB |
| --- | ---: | ---: | ---: |
| Force along initial hybrid Phe1 CA direction | 3.341027 | 51.063006 | 0.798547 |
| Force along full deformation chord | 3.638430 | 19.928221 | 2.466204 |
| Maximum atom force | 8.204600 | 51.063006 | 2.140905 |
| RMS atom-vector force | 3.666509 | 9.121235 | 0.921311 |
| RMS force difference from reference | — | 10.141916 | 3.511987 |

This initial force comparison cannot determine the energy cost of the complete
motion. It does not identify a unique MM parameter defect or validate MACE.
The remaining three cells retain the original fixed geometries, method and
complete four-cell denominator. Root's measured-cost decision is recorded in
`workspaces/lanm_coupled_mechanics_20261017/lbt_whole_reference_continuation_v1/`.
Post-scout execution review and actual saved-data validation precede those runs.

## Native evidence and limits

ORCA terminated normally after18 reported SCF cycles and produced the full
gradient, including ECP and D4 derivative evidence. Charge−2, La effective
singlet,1030 explicit electrons, vacuum scope and12MPI/MaxCore3100 match the
declared recipe. Final energy is−6765.679869776899Eh.

Orbital-gradient and rotation criteria pass. Printed energy-change, maximum
density-change, RMS density-change and DIIS-error flags remain false under
the original native convergence semantics. No new all-flags acceptance rule
was imposed. Native analytic-gradient admission is not an independent numerical
derivative qualification or a proof of electronic stability.

Conditional free termini, the vacuum dianion, post-added solvent, frozen-f/SOC
limitations and nonstationary endpoints remain explicit in every result.
No finite deformation work, La/Dy preference or production promotion is established.

## Provenance and cost

- Scout13051:23,741 scheduler seconds on12CPUs,284,892 allocated CPU-seconds.
- Saved collector13283:3seconds on1CPU,3allocated CPU-seconds; zero molecular calls.
- Exactly one new DFT call; ALPB and both model comparisons were reused.
- `SCOUT_RESULT.json`, `COLLECTION_RESULT.json`, `FORCE_COMPARISON.json`,
  `NATIVE_RECEIPT.json`, `SACCT.txt` and `ROOT_REVIEW.json` preserve actual evidence.
- Independent reviewer checked raw pins, native termination,256 gradients,
  source admission rules and the interpretation; no declared-admission blocker.

Three more runs at the scout's cost would allocate854,676CPU-seconds. This is
a planning comparison, not a runtime forecast, imposed limit or measured cost.
