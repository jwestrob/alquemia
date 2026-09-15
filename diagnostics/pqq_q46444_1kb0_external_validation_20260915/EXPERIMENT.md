# Q46444 / 1KB0 fixed-core-v3 external validation

Date preregistered: 2026-09-15  
Protocol: `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`  
Target: Q46444 quinoprotein alcohol dehydrogenase crystal structure 1KB0

## Question and frozen answer rule

This is a one-target external structural-transfer test of the already calibrated
PQQ-MDH La/Ca discriminator. It asks whether an independently selected,
Ca-containing quinoprotein site lands in the released Ca-supported score band.
It does not refit chemistry, atom selection, the aquo gauge, or either band.

The expected result was frozen before ORCA execution:

- `S <= 14.857129202922806 kcal/mol`: **external validation pass**
  (`Ca-supported`).
- `14.857129202922806 < S < 23.460061205609236 kcal/mol`:
  **indeterminate and not a pass**.
- `S >= 23.460061205609236 kcal/mol`, an invalid arm, or a paired-arm
  invariant failure: **external validation failure**.

Here `S = [(E_Ca - E_La) - delta_E_aquo] * 627.509474`. The locked aquo
reference and exact released bands come only from the immutable calibration
result. No result from 1KB0 may alter them.

## Exact source and selectors

- Raw source: `/groups/banfield/projects/environmental/sr/srvp2020/protenix/mmcif/1kb0.cif`
- Raw-source SHA-256:
  `2b26ccaecd3ed6c89af5e6d8404e39ee09787d79d1904c7c99c48658685aa916`
- Resolution: 1.44 A; model 1; chain A only; no selector fallback.
- Native metal: `A:CA801/CA`, occupancy 1.00, B = 21.65 A2, coordinate
  `(19.844, 56.738, 15.434) A`.
- PQQ: `A:PQQ1800`, complete canonical CCD PQQ, frozen as PQQ(3-).
- Protein roles: anchor Glu `A:GLU185`; anchor Asn `A:ASN263`;
  catalytic Asp `A:ASP308`; D+2 homolog `A:THR310` (recorded but excluded);
  second-shell cationic partner `A:LYS335`.

The frozen typed coordination number at 3.1 A is 7: PQQ O7A/N6/O5,
GLU185 OE1/OE2, ASN263 OD1, and ASP308 OD1. Fragment membership is
role-based and must not change if a distance is recomputed after protonation.

## Dry-selection policy

Preparation starts from the exact hashed raw crystal and retains model 1,
chain A standard amino acids, the selected Ca coordinate, and the selected
PQQ only. It removes all 1,016 source waters and every noncore heterogen,
including TFB1810, HEC802, TRO512, and all glycerols. No synthetic water,
counterion, buffer molecule, or replacement ligand is added. In particular,
`A:TFB1810/OXT` is 3.338 A from Ca: it lies outside the frozen 3.1 A donor
gate but inside 3.6 A, and is deliberately removed without replacement.

The native Ca atom name, element, and coordinate are retained through the
standard-residue-only protonation step; only its residue label is normalized
to `LA` so the immutable v3 machinery can process it. The final La and Ca QM
arms use byte-identical nonmetal coordinates and differ only in metal identity,
total charge, and electron count. There is no geometry relaxation, point-charge
embedding, or source/synthetic water in this experiment.

## Interpretation boundary

1KB0 is a Ca-containing QH-ADH structure and an independent structural target
relative to the frozen 25-member calibration. Its label is biochemical and
structural support for a Ca-class test, not a direct matched La/Ca equilibrium
constant. Passing therefore supports transfer of the discriminator to this
site; it does not by itself estimate affinity or establish physiological metal
exclusion.

