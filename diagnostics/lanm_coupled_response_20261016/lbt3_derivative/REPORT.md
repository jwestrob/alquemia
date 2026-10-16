# Solvent-containing forces pass the real donor-direction check

The declared intact Asp7 carboxylate rotation in La-derived LBT3 passes the
unchanged numerical criterion for both metals and both displacement sizes.
All 8 MACE and 16 standalone GFN2 fixed energies were admitted. Coordinates,
proton inventory, chemical state and component definitions were unchanged.

| Metal | Analytic derivative | Numerical at 0.001 Å | Numerical at 0.0005 Å |
|---|---:|---:|---:|
| La | −8.16533177 | −8.16532250 | −8.16532412 |
| Dy | +7.02280586 | +7.02280463 | +7.02280856 |

Units are kcal/mol/Å of maximum physical chord displacement. The largest
composite residual, including the paired derivative, is about 0.0000105.
Individual MACE, GFN2-vacuum, GFN2-water and solvent residuals and both-step
comparisons remain in RESULT.json. The fixed tolerance was 0.2 plus 5% of the
analytic magnitude; no threshold was changed after viewing these results.

The opposite slopes confirm a metal-dependent local response in this model.
This tests numerical consistency along one real motion, not full force-field
accuracy, Dy electronic-state accuracy, affinity or protein-level transfer.

## Executed work

Force/geometry staging12512 used14s/oneCPU. Finite MACE12519 used10s/12CPUs
and oneH200. Finite solvent12520 used100s/fourCPUs; measured totalCPU367.902s.
Saved comparison12521 used1s/oneCPU. These four allocations add535CPU-seconds
and10GPU-seconds. Branch total through12521 is1321CPU-seconds/48GPU-seconds,
including previously reported preparation and failures; not measured utilization.

The passed check permits the separately frozen LBT3_ACCOMMODATION_PLAN.md.
No optimization has yet run. Root owns submission of the first7CCO/La scout,
followed by the unchanged other three searches if execution/physical admission
passes, independently of the preferred sign or amount of energy reduction.
