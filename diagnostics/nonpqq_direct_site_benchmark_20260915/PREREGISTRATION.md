# Frozen non-PQQ direct-site La/Ca benchmark

Date frozen: 2026-09-15
Status at freeze: input preparation only; no ORCA calculation submitted or run
Panel manifest: `panel_manifest.json`

## Question

Can the architecture-neutral, frozen-site La/Ca vertical-exchange score carry
information outside the PQQ dehydrogenase motif/core-charge split?

This panel deliberately contains six crystallographically explicit, CN7 Ca
sites from three unrelated proteins. It is **not** used to refit the PQQ
protocol or its 19.1586-kcal/mol threshold. The first calculation made after
this freeze is confirmatory with respect to the choices below; no donor,
water, site, or protein may be dropped because its score is inconvenient.

## Frozen structural protocol

1. Use the first deposited model of the exact source mmCIF named and hashed in
   `panel_manifest.json`.
2. Add standard-residue and crystallographic-water hydrogens at pH 7 with
   `pdbfixer_standard_residue_protonation_v2`. Do not build missing residues.
   Source heavy atoms are retained; their maximum coordinate displacement is
   measured in the preparation report.
3. Select each Ca atom by the exact author chain, residue number, insertion
   code, and atom name in the panel manifest. Automatic site ranking is
   forbidden.
4. Apply `generic_vertical_exchange_native_r2scan3c_v2`: typed CN is measured
   at 3.1 A and chemically complete donor fragments are included through
   3.3 A. Every listed first-shell water is included as complete neutral H2O.
5. Generate a vertical La(III)/Ca(II) pair at the **same deposited Ca
   coordinate**. All nonmetal atoms and coordinates must be byte-identical
   between the two XYZ files. There is no geometry relaxation and no
   metal-specific donor or water selection.
6. Electronic protocol, if executed later: ORCA 6.1.1, native r2SCAN-3c,
   `NoAutostart CPCM(Water) DefGrid3`, native composite basis/ECP policy, and
   symmetric CN8 aquo reference
   `aquo_cn8_symmetric_vertical_native_r2scan3c_cpcm_v2`.

The score is

```text
S = [(E_Ca,site - E_La,site) - (E_Ca,aquo - E_La,aquo)] * 627.509474
```

Positive means La-favoring relative to the fixed aquo gauge; negative means
Ca-favoring. This is a frozen electronic compatibility score, not a binding
free energy or a predicted Kd.

## Frozen evidence and interpretation

### GGR/MglB (`P0AEE5`, 1GLG)

This is the sole primary pass/fail control in this six-site panel. Ca and La
were compared in the same assay: Kd(Ca) = 25 +/- 11 uM and Kd(La) = 729 +/- 4
uM, about 29-fold Ca preference (DOI `10.1021/bi00468a021`). The preregistered
direction is therefore `S < 0` for the single A:CA312 site. Failure cannot be
rescued by any aequorin or parvalbumin result.

### Aequorin (`P07164`, 1SL8)

The reported apparent association constants, 1.92 and 1.38 uM^-1 for Ca and
La, respectively, are a weak **protein-level** Ca preference (PMID 2153542).
The measurement cannot be assigned honestly to one of the three deposited EF
hands. The primary output is consequently the ordered vector
`[EF1, EF3, EF4]`, never a selected best site. The median, range, and sign count
are frozen descriptive summaries only. Aequorin has no scalar calibration
label and cannot by itself pass or fail the discriminator.

### Carp parvalbumin pI 4.25 (`P02618`, 4CPV)

The La study assigns Kd(La) = 20 pM to the CD site and 48 pM to the EF site
(DOI `10.1021/bi00294a030`). Published Ca affinity comes from a different
experimental context, so neither site is promoted to a same-assay quantitative
truth label. The primary output is the ordered vector `[CD, EF]`; both scores
are reported separately. These are supporting La-positive observations, not
threshold-setting controls, and their sub-kcal ordering is not a valid target.

## Water and geometry limitations fixed before scoring

1GLG is a dry protein-CN7 site. Each 1SL8 site and the 4CPV EF site contains
one crystallographic first-shell water; the 4CPV CD site is protein-CN7 and
dry. X-ray structures resolve water oxygen but not hydrogen orientation. The
primary panel uses exactly one deterministic PDBFixer orientation and records
its coordinates and hash. Water-orientation sensitivity may be run only as a
separately labeled follow-up and may not replace the primary result.

All six states are conditioned on deposited Ca geometry. A negative result can
therefore expose failure of the frozen electronic/site model, but a positive
result does not establish that reorganization, cooperativity, entropy, proton
linkage, or absolute binding affinity has been reproduced.

## Falsification and prohibited reinterpretations

- A nonnegative GGR score is a failure of the primary non-PQQ direction test.
- The three aequorin sites must remain an ordered set; selecting whichever site
  agrees with the weak global assay is prohibited.
- The two parvalbumin sites must remain separate; their cross-study labels may
  not override a GGR failure.
- No score from this panel sets or validates a universal numerical threshold.
- No result here demonstrates within-lanthanide discrimination.
- A preparation with missing/extra typed donors, a missing listed water,
  moved source heavy atoms, or non-identical La/Ca nonmetal coordinates is
  ineligible and must not be submitted.
