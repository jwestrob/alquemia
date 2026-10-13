# All 15 two-Ca sources are prepared for matched scoring

The spectator-aware preparation passed for **15/15 sources**, with no exclusions.
Both native Ca ions and every source heavy atom remain present at their original
coordinates. The new local contexts match the original canonical chemical state
and physical-map requirements; no donor, proton, water, cutoff or band changed.

| Protein | Prepared | Full-source heavy atoms | Paired context atoms | Spectator distance to nearest retained heavy atom |
|---|---:|---:|---:|---:|
| Q9Z4J7 | 5/5 | 4,849 | 192 | 23.79–23.99 Å |
| C5AXV8 | 5/5 | 4,543 | 192 | 24.02–24.21 Å |
| Q88JH5 | 5/5 | 4,922 | 192 | 23.74–23.81 Å |

All source-heavy coordinate residuals are exactly zero. The active-site selector
is taken from each real audit: Q9 sample3, C5 sample1 and all Q88 samples use
chain D's PQQ Ca; the other sources use chain B. The other ion remains the
spectator. All spectator distances exceed the frozen 4.3 Å retained-heavy gate,
and none appears in local atom/cap/bond maps. All paired contexts preserve their
original canonical state/protein signatures and source mappings.

The new narrow protonation wrapper also passes an actual old one-Ca regression:
P12293 sample0 reproduces archived core/context XYZ and expanded maps/state
exactly, with no displacement and identical XYZ hashes. The remaining ten
preparations reused that passed receipt; they did not repeat the regression.

## What this establishes

The original one-Ca failures were not a reason to omit the structural Ca silently.
The explicit two-Ca sources now pass the same local chemistry and mapping gates,
so their effect can be tested fairly with the existing scorer. These are distinct
occupancy hypotheses derived from actual refolds, not repaired old observations
or independent biological replicates. The original 125-source benchmark remains
unchanged.

This stage establishes preparation coverage, **not improved classification**.
The spectator remains in the full source but is excluded explicitly from the
local Hamiltonian; its global electrostatic/mechanical contribution is not
modeled. Separate native/common-pool and corrected-solvent scoring must determine
whether the new structures actually improve discrimination under frozen bands.
No energy endpoint or model forward was run by this preparation adapter.

## Actual execution

- Job2969: five Q9 preparations plus one-Ca regression, normal completion,
  28 scheduler seconds × 6 CPUs; preparation wall 26.54 s, process CPU 81.25 s,
  batch MaxRSS 738,908 KiB.
- Job2985: ten C5/Q88 preparations, normal completion, 30 scheduler seconds ×
  6 CPUs; preparation wall 28.99 s, process CPU 136.79 s, batch MaxRSS 633,520 KiB.
- Combined allocation: 348 CPU-seconds. Zero model/endpoint calls. All scientific
  preparation and checks ran in Slurm; login work was metadata, code and syntax.

Exact result pins, all 15 rows, preservation/boundary receipts and measured
costs are in `FINAL_RESULT.json`. Authoritative downstream artifacts:

- `workspaces/ca_coverage_recovery_20261012/spectator_preparation_q9_v2/execution/FINAL_PREPARATION.json`
- `workspaces/ca_coverage_recovery_20261012/spectator_preparation_c5q88_v1/execution/FINAL_PREPARATION.json`

The code and consequential settings are identical across both jobs. Original
helpers and production remain unchanged. Root owns downstream scoring and any
promotion decision.
