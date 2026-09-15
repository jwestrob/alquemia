# Preregistered fixed-core PQQ-MDH calibration

**Protocol ID:** `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`

**Status at freeze (2026-09-14):** design only. No v3 electronic energy has
been calculated or inspected. This document fixes the chemistry, calibration
rule, and pass/fail criteria before input preparation and execution.

## Question

Can one chemically explicit, fixed PQQ-MDH active-site core distinguish the
11 previously verified lanthanide-dependent PQQ-MDHs from the 14 previously
verified calcium-dependent PQQ-MDHs, without changing the carve between
proteins or choosing a cutoff after seeing the new energies?

The experiment deliberately recalibrates before testing matched crystal
structures. The original 25-member panel is calibration data, not new external
validation. PDB 1H4I (MxaF) and PDB 4MAE (XoxF) are reserved as the primary
structure-level transfer test and must not be calculated under v3 unless the
calibration gate passes.

## Frozen calibration panel

Use every entry in `testset_expansion/frozen_v0_ids.tsv`: 11 Ln-labeled and
14 Ca-labeled proteins. No protein may be removed, relabeled, replaced, or
weighted. All 25 must produce one valid La/Ca pair.

Use the preserved source CIF for each protein, normalize metal/PQQ identities,
and protonate by the repository's current deterministic PQQ workflow. Retain
the source geometry; do not refold, optimize, minimize, add a water, or select
an alternative conformer. La and Ca calculations for a protein must have
byte-identical nonmetal coordinates and differ only in metal identity, charge,
and electron count. Both are closed-shell singlets.

The source structures contain no resolved waters. The primary v3 state is
therefore explicitly dry; synthetic waters are forbidden.

## Frozen QM core

The QM coordinates contain full PQQ(3-), the metal at the source metal
coordinate, and side-chain fragments from these complete, residue-identity
selectors on chain A:

1. the conserved anchor Glu;
2. the conserved anchor Asn;
3. the catalytic Asp immediately following the conserved Trp;
4. the residue two positions after that catalytic Asp **only when it is an
   acidic direct donor** (Asp in all 11 Ln controls; Ala/Ser/Thr in all 14 Ca
   controls are recorded but not carved); and
5. the conserved cationic partner hydrogen-bonded to the catalytic Asp (Arg in
   23 proteins, Lys in 2).

The fifth fragment is a required second-shell electrostatic partner, not a
metal donor: its closest N is 4.329--5.102 A from the metal but 2.584--3.333 A
from catalytic Asp OD2. It must never inflate coordination number.

Selectors, source hashes, fragment atom identities, link hydrogens, formal
charges, and final coordinate hashes must be recorded before execution.
Distance-radius membership is not allowed to alter this core. The expected
total La/Ca charges are -2/-3 for Ln controls and -1/-2 for Ca controls.
Any mismatch is a preparation failure.

For descriptive QC, geometric CN is the number of PQQ/protein O or N atoms
within 3.1 A of the metal. Fixed-core inclusion does not itself count as
coordination. Calibration admission requires source CN >= 6; the production
prospecting gate remains CN >= 7 and is not changed by this experiment.

## Frozen electronic model

- ORCA 6.1.1 single points.
- `r2SCAN-3c NoAutostart CPCM(Water) DefGrid3`.
- Native r2SCAN-3c basis/ECP policy for both metals; no custom La override.
- No point-charge embedding, MM environment, geometry relaxation, BSSE
  correction, or post hoc energy correction.
- Use unrounded final single-point energies; any abnormal termination,
  nonconvergence, charge mismatch, coordinate mismatch, or missing selector is
  invalid and may only be rerun unchanged.

The finalized symmetric CN8 aquo record is a reporting gauge only:
`reference_inputs/aquo_cn8_native_r2scan3c_v2/aquo_reference.json`, frozen
SHA-256
`cf4f42ca1c768ed55dedbb2360f21eb1106ca77bada4347e8db9e3511836f0d2`.
It must not be edited to declare v3 compatibility. Any future compatibility
extension requires a new immutable reference record.

## Frozen score and calibration rule

For protein `i`, define the primary raw site contrast

`R_i = (E_Ca,site - E_La,site) * 627.509474 kcal mol-1`.

Larger `R` is more lanthanide-like. If `A` is the analogous aquo contrast,
the familiar exchange score is `S_i = R_i - A`; this subtracts the same
constant from every protein. All classification margins and conclusions are
therefore made in `R`, where aquo compatibility cannot affect the answer.

After all 25 valid scores exist, define

- `U = max(R_i)` among the 14 Ca controls;
- `L = min(R_i)` among the 11 Ln controls; and
- only if `L > U`, `T_R = (U + L) / 2`.

The protocol is **CALIBRATABLE** only if all four conditions hold:

1. 25/25 proteins have valid, unchanged v3 La/Ca pairs;
2. `L > U` (complete directional separation; AUROC = 1);
3. the empty-class gap `G = L - U` is at least 5.0 kcal/mol; and
4. leave-one-protein-out midpoint classification is correct for 25/25 proteins,
   with each omitted protein classified using the midpoint derived from the
   other 24.

The 5.0 kcal/mol floor was fixed before v3 energies from the pre-existing
maximum cross-preparation movement of approximately 4.49 kcal/mol, rounded
upward. If any gate fails, v3 is **NOT CALIBRATABLE**: do not optimize a Youden
cutoff, remove an outlier, shift the aquo gauge, or consume the reserved
crystal pair as a rescue. A chemistry revision requires a new protocol ID and
a complete 25-protein rerun.

If the gate passes, release the midpoint threshold and the conservative bands:
Ca-supported at `R <= U`, Ln-supported at `R >= L`, and indeterminate in the
unoccupied interval `U < R < L`.

## Reserved structure-level transfer test

Only after the calibration gate and threshold are locked, prepare and run:

- PDB 1H4I MxaF as the Ca-dependent negative control; and
- PDB 4MAE XoxF as the Ln-dependent positive control.

Both must use exactly the frozen v3 selectors and electronic model. The
transfer test passes only if both are valid, 1H4I has `R <= U`, and 4MAE has
`R >= L`. A score in the calibration gap is indeterminate; a wrong-band or
unscorable result fails. PDB 6OC6 may be reported only as a secondary geometry
check because its sequence is already represented by C5B120 in calibration.

The crystal test is a held-out structural transfer test, not a blind or
sequence-independent validation: mature 1H4I occurs within the P16027
calibration sequence, and prior protocol versions of these structures are
known. All such limitations must remain explicit in the result.
