# Second-shell donor context and PQQ hydrogen preparation

**Both pilots completed. Explicit context modestly increases the already correct
alpha/GGR separation; neither pilot improves the PQQ gap or demonstrates new
classification accuracy.** Production, released references and the default
scorer remain unchanged.

## Fixed-coordinate context

The same geometric rule adds complete polar peptide/side-chain units within
3.5 Å of direct donor functional groups. Source connectivity merges overlapping
fragments. PQQ, shared source coordinates, surviving caps, water identity and
water inventory remain exact. Alpha starts with the previously improved water
hydrogens in both representations. No geometry search or new electronic method
is mixed into this comparison.

| Expected La-like minus Ca-like | DFT core | DFT context | Native MACE core | Native MACE context |
|---|---:|---:|---:|---:|
| XoxF 4MAE − MxaF 1H4I | 30.572 | 28.179 | 102.222 | 89.825 |
| alpha 1F6S − GGR 1GLG | 10.924 | 14.815 | 8.304 | 35.437 |
| alpha 6IP9 − GGR 1GLG | 14.573 | 18.422 | 9.613 | 28.417 |

Units are kcal/mol; positive is the expected relative direction. The two alpha
structures are replicas of one biological comparison. PQQ functional association
and condition-qualified alpha/GGR affinity evidence remain separate. Every case
was already inspected during earlier development. No threshold was fitted or
inherited, so these are relative contrasts, not new absolute classifications.

Native MACE agrees on the direction of all three *gap changes*, but exaggerates
them. Its individual same-site changes disagree in sign with CPCM-DFT for GGR
and both alpha structures. It also increases variation between the alpha
replicas. Native OMOL is a vacuum model with total-charge conditioning; the DFT
target uses CPCM. The result does not justify adding its context difference to
CPCM-DFT as a quantitatively validated environmental correction, or rehabilitate
previous failures of native whole-protein charge conditioning.

Expanded sizes are 154, 202, 125, 106 and 104 atoms, respectively. Added formal
charge is zero in both PQQ contexts and −1 in GGR and both alpha contexts. Thus
the difference in total charge between the two members of each biological
comparison does not change. Composition, caps and CPCM cavity still change;
the contrast shifts cannot be uniquely attributed to hydrogen bonding.
This is not a partition-invariance test on a common fixed physical boundary.

All 10 new native DFT endpoints terminated normally. Four MACE calls per structure
supply the exact original/expanded Ca/La comparison: 20 completed calls. DFT
cost 2,244 seconds on 64 CPUs; MACE cost 145 seconds on one GPU and 16 CPUs.
The larger explicit DFT context is useful for this development test but expensive
for a routine scanner. Cheap MACE execution alone does not establish useful
predictive accuracy.

Unrounded energies, matched shifts and component diagnostics are in
[result_v1/result.json](result_v1/result.json); [tables](result_v1/TABLES.md).
See [the preserved collection field-name note](COLLECTION_NOTE.md).
SCF contains the CPCM contribution: the separately reported CPCM dielectric term
is an overlapping diagnostic, not another additive energy.

## PQQ physical-H follow-on

Four native MACE searches move physical hydrogens in the same PQQ contexts,
with heavy atoms and artificial caps fixed. Transferring only original-core
hydrogens back to the original 73/80-atom cores lowers all four native DFT endpoint
energies by 170–263 kcal/mol. The Ca/La difference mostly cancels: the PQQ gap
changes 30.572 → 30.072 kcal/mol. Native MACE instead changes 102.222 → 104.027.
This establishes retained ordering and lower electronic energies, not improved
DFT discrimination. See the [separate report](PQQ_H_REPORT.md).

The frozen 0.35 Å domain is active for 16–35 hydrogens per endpoint. All four
meet the declared constrained force criterion; one hits the 200-iteration solver
limit. No unconstrained minimum, entropy correction or new calibrated band is
claimed. The same final-proposal rule was used regardless of the result.

## Delivery and interpretation

- New static protocol: `native_r2scan3c_cpcm_fixed_donor_second_shell_v1`;
  matched native MACE protocol: `native_OMOL_second_shell_fixed_v1`.
- New preparation/scoring protocols:
  `PQQ_native_OMOL_context_physical_H_preparation_v1` and
  `native_r2scan3c_PQQ_context_prepared_physical_H_v1`.
- Eleven real-fixture checks pass with zero skips. These include actual MACE/DFT
  receipts, exact source/PQQ/unchanged-H checks, graph overlap and charge parity,
  corrupted-neighbor rejection, constrained hydrogen mapping and score algebra.
- Total new work: **14 DFT single points**, 20 static MACE endpoint calls,
  four MACE hydrogen searches with 2,113 recorded objective invocations, and
  four proposed-core MACE calls. ASE may cache repeated coordinates.
- Total allocation: **171,088 core-seconds and 409 GPU-seconds** across four
  completed jobs. Login preparation/tests and reused earlier calculations are
  additional. Exact scheduler accounting: [COSTS.json](COSTS.json).
- A first preparation failed JSON integer-key replay checking before execution.
  Its artifacts remain preserved; the fix changed no scientific input.

**Recommendation:** retain the current PQQ baseline and prioritize the separately
successful cheap water-preparation scanner track. These two second-shell pilots
supply useful constraints on the mechanism, but do not earn a PQQ classifier
promotion. Larger canonical evaluation remains possible under a separately
frozen preparation policy; no such expansion was launched here.

[Approved static scope](AGREEMENT.md), [approved H scope](PQQ_H_AGREEMENT.md),
[commands](COMMANDS.md). The vault note is
`agent-captures/2026-09-19_laca-second-shell-fixed-context.md`.
