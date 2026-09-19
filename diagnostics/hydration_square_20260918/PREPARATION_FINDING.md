# Stretched water hydrogens in both alpha-lactalbumin parents

Read-only inspection completed before new electronic calculations. The source
cores, archived energies and production preparation are unchanged.

| Source | Water | O–H1 (A) | O–H2 (A) | H–O–H (degrees) |
|---|---|---:|---:|---:|
| 1F6S | A:HOH211 | 1.182354008 | 1.162862417 | 105.705557 |
| 1F6S | A:HOH212 | 1.173068625 | 1.176653730 | 105.488501 |
| 6IP9 | A:HOH310 | 1.185938447 | 1.182476638 | 105.873357 |
| 6IP9 | A:HOH322 | 1.189716773 | 1.180331309 | 105.343009 |
| 6IP9 | A:HOH326 | 1.172886184 | 1.187677566 | 105.514902 |

The real source atom mappings, input/output hashes and archived endpoint checks
are retained in `workspaces/hydration_square_20260918/source_audit_v2.json`.
Both source protonation manifests name PDBFixer 1.12.0 / OpenMM 8.5.1, pH 7.
`protonate_cif.protonate` calls `fixer.addMissingHydrogens(ph)` without a supplied
forcefield. Inspection of the installed implementation confirms that PDBFixer
forwards `forcefield=None` to `Modeller.addHydrogens`. Its fallback uses generic
hydrogen springs, a repulsive nonbonded potential and a water angle term, then
minimizes the added hydrogens. This does not constitute a physical water model
or quantum relaxation of the pocket. The observed distortion is already present
in the archived protonated structures, before carving.

## Proposed contained repair

Use the geometry of the existing pinned three-atom `bulk_water.xyz` as a common
internal-water template (O–H approximately 0.9572 A, angle 104.52 degrees).
For every retained water in both parents, keep O fixed and preserve the original
water plane and bisector; change only its two hydrogen positions. Keep protein,
caps, protonation, charges and all other coordinates byte-identical. This is a
deterministic internal-geometry normalization, not a physically optimized water
orientation or proof of preferred occupancy. It introduces no metal-dependent
coordinate rule.

Recompute both metals on each repaired full-water parent (four endpoints), then
remove each individual water in turn (ten endpoints). Old full-water endpoints
separately quantify the normalization's effect. All five water-removal squares
are reported; none is selected from its score. The repair's effect and hydration
effect must remain separate in the report.

Jacob's geometry choice is pending. No new DFT or MACE evaluations, allocations
or classifier fits have run. Four real-fixture software tests pass in 0.341 s:
source water identities/neutral deletion, preservation of heavy atoms/planes/
bisectors and template geometry, explicitly corrupted mapping rejection, and
actual archived output/receipt/component parsing. These tests do not validate
water thermodynamics.
