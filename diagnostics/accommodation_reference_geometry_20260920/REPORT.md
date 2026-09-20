# Known-reference fold donor geometry and mapping support

**62 of 233 prepared reference folds trigger the existing <2.2 Å contact screen;
all 62 support the current physical terminal-carboxylate mapping.** They span
14 protein groups, not 62 independent biological observations. These belong
to the consumed PQQ functional-reference set, not new blind affinity labels. No energies,
classification labels or scores were used to select them. No molecular
calculation, displaced endpoint, ion movement or source-coordinate change ran.

| Source conditioning | All folds | Prepared | Unsupported | Flagged Glu | Flagged extra Asp | Any flag |
|---|---:|---:|---:|---:|---:|---:|
| Ca | 125 | 110 | 15 | 19 | 1 | 20 |
| La | 125 | 123 | 2 | 42 | 0 | 42 |

Every one of the 250 original statuses is retained in the detailed inventory;
the 17 unsupported preparations have unavailable screen results, not passes.
The predicate is strictly below 2.2 Å, reused from the existing PLM geometry
review. It is a La-oriented warning, not a universal bond-validity cutoff; a
Ca-conditioned structure is not defective merely because this screen fires.

## What the test set can address

All 42 La-conditioned flags concern the anchor Glu (nearest O 2.0520–2.1939 Å).
None has an extra-Asp oxygen below the cutoff. Among Ca-conditioned structures,
the 19 Glu flags span 2.1178–2.1974 Å; the one extra-Asp flag is
`mmol_1770-pqq-la_model__conditioned_Ca__seed-1_sample-1` at 2.1016 Å.
These contacts are substantially less compressed than the two PLM extra-Asp
origins (1.766/1.847 Å).

This supplies known-label structural cases for a possible **Glu-response** test,
with source-metal conditioning retained. It supplies very little replication of
the exact **extra-Asp compression** phenotype motivating the ongoing PLM pilot.
The pending extra-Asp DFT checks cannot automatically qualify Glu accommodation.
A supported coordinate map establishes that the intended motion is implementable;
it does not establish a good energy surface, better classification or a relaxed
structure. No follow-on energy manifest or case selection was launched here.

## Source and implementation checks

Source is the root's immutable `folds_v1/preparation_reconciled_v1.json`.
Metal coordinates agree with its recorded actual source ion. Ca/La prepared
coordinates are paired exactly, apart from metal identity; charges differ by one.
Each oxygen maps back to the corresponding original fragment atom record.
Non-Asp residues at the extra-acidic homolog position are explicitly
not applicable; absent nonacidic fragments are not invented.

For each flagged structure, the existing source graph and complete-context
preparation reconstruct the physical mapping at its unchanged origin. Terminal
Glu CG–CD or Asp CB–CG modes retain both carboxylate oxygens, with existing bond,
cap, Jacobian and fixed-scaffold checks. The same Cartesian coordinate measure
applies to both metal endpoints. All 62 mappings pass. Internal algebraic mapping
checks create no scientific endpoint or altered source file.

## Artifacts and cost

- [SUMMARY.json](SUMMARY.json): compact counts, ranges and final source pins.
- [INVENTORY.json](INVENTORY.json): every source condition/status, role and distance.
- `workspaces/accommodation_reference_geometry_20260920/mappings_v1/result.json`:
  all mapping outcomes; source graph mappings reside in its `maps/` directory.
- `scripts/accommodation_reference_geometry.py`: runnable `scan` / `mappings`.
- Job1203826 performed mapping checks only, in 2.328 process seconds using the
  existing CPU environment; [COSTS_sacct.tsv](COSTS_sacct.tsv) records actual
  allocation: 2 whole Slurm wall seconds ×64 CPUs =128 allocated core-seconds
  (whole-second resolution); recorded CPU time55.590 s. Zero reported batch
  MaxRSS is unavailable, not a measured zero. No GPUs, ORCA or MACE calls. Local inventory/summary work is outside
  that allocation receipt.

The fixed cutoff and selection use no label or score. Production remains unchanged.
