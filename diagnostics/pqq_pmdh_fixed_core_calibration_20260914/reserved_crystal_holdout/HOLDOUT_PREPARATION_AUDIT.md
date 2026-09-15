# Reserved PQQ-MDH holdout preparation audit

Date: 2026-09-15
Protocol: `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`
Implementation commit: `df166c2`

## Release decision

**READY FOR ORCA:** the atomic primary pair 1H4I + 4MAE prepared as four
La/Ca single-point legs. The immutable calibration gate was independently
rechecked as `CALIBRATABLE` (25/25 valid, AUROC 1.0, LOO 25/25, frozen S bands
Ca <= 14.857129203 and Ln >= 23.460061206 kcal/mol).

6OC6 was not released. Its attempted standard-AA protonation would reconstruct
130 unresolved heavy atoms across 44 noncore residues; both primary structures
require zero repairs. The zero-heavy-repair gate was retained, and the
nonindependent secondary was recorded as `secondary-not-run`. The invalid,
never-executed attempt is preserved at:

`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/preparation_failed_secondary_missing_heavy_atoms_20260915T050532Z`

## Audited preparation

Preparation manifest:

`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/prepared/holdout_preparation.json`

SHA-256: `c7942ce190b1af50617798aa2ac6ca2aee40e26602c54b34217ee2885375470f`

| PDB | Biological role | Typed CN | O/N | Fragments | QM atoms | La/Ca charge | Heavy repairs | Max heavy displacement |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 1H4I | Ca/MxaF primary | 6 | 5/1 | 5 | 72 | -1/-2 | 0 | 0.000 A |
| 4MAE | Ln/XoxF primary | 9 | 8/1 | 6 | 79 | -2/-3 | 0 | 0.000 A |

For both targets, retained source heavy-atom identities and coordinates are
exactly preserved, La/Ca arms have byte-identical nonmetal coordinates, and
the arms differ only in metal identity, total charge, and electron count.
Each has zero unexpected typed donors, nearby untyped O/N atoms, or nearby
sulfur atoms under the frozen checks. PQQ is complete `pqq_ox_3minus_v1`;
there are no source or synthetic QM waters, no point charges, and no geometry
relaxation. The 4MAE coordinating `A:15P603/OXT` crystallization adduct is
explicitly excluded without replacement, as preregistered.

No ORCA output, energy, or holdout score existed when this audit was written.
