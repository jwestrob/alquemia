# Uniform strict native scalar qualification

Approved by root under Jacob's overnight discretionary authorization, 2026-09-23. The latest coordination message supersedes the unexecuted cold-seed/pass2-seed full-panel design and separate ten-cell direct-start check. Neither superseded design was prepared or submitted.

## Question and population

Does native GFN2 with TolE=1e-10 Eh give repeatable scalar energies from a fresh atomic guess and an archived loose-cold seed, while retaining useful fixed-pool PQQ discrimination? Use the exact 32 sources/384 cells from `workspaces/precision_pool_continuation_20260923/INVENTORY_v2.json`: 25 canonical calibration references, three consumed crystals, and four consumed noncanonical folds. Same original three geometries (origin, adaptive_Ca, adaptive_La), two metals and vacuum/ALPB water. No new geometry, MACE, DFT, optimization or labels.

## Fixed two branches

A: 384 fresh native singlepoints, explicit NoAutostart, no seed files.
B: 384 strict continuations from each cell's exact original loose-cold GBW+xtbw pair. Reuse the ten exact branch-B outputs from strict20; 374 new B calls. Total 758 new / 10 reused evaluations. The earlier pass2-seeded strict20 results remain supplementary and never substitute for A or B.

Both: native GFN2, native mixer, electronic temperature 300 K, MaxIter500, TolE1e-10 Eh, original solvent defaults/parameters/charge/multiplicity/coordinates, one rank. All new tasks remain on the existing finite ORCA manifest runner. Two separate contained manifests/locks execute sequentially in one 32-worker, 32-CPU, 64-GiB CPU-only allocation on the approved GPU-partition host. No retries, third starts or outcome-selected stage.

## Frozen gates and reporting

Record every cell, failure and actual guess; confirm effective TolE/native mixer/state/parameter export and actual XTBRESTART only for B. Require |A-B| <=0.1 kcal/mol per cell; same-geometry composite Ca−La contrasts and each pooled mathematical/operational contrast must agree within 0.2 kcal/mol. Operational candidate selection retains the existing 0.1 kcal/mol origin preference rule. Both branches keep raw pools, candidate choices, components, old-band transfer and their own calibration records. Qualified references require all 25 designated calibration sources and all numerical gates, using the unchanged extrema/minimum-gap rule. Crystals and four folds never fit bands. Missing values remain null, never another branch's output.

A is the proposed one-call method and must succeed independently. B uses a prior loose-cold call plus this continuation; report both cost components, never call it one-call production. Development validation computes both branches. Energy repeatability does not qualify analytic forces, equilibrium populations, a unique electronic solution or broad biological accuracy. No production/default change.
