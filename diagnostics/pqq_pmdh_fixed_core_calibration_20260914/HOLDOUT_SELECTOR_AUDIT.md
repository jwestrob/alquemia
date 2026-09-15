# Reserved PQQ-MDH crystal holdout selector audit

Date: 2026-09-14  
Protocol: `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`

The exact machine-readable selectors and source hashes are frozen in
`holdout_spec.tsv`. This audit only examined coordinate metadata. It prepared
no QM input and inspected or calculated no v3 energy.

## Frozen disposition

- **Primary negative:** 1H4I, Ca-dependent MxaF, source chain A.
- **Primary positive:** 4MAE, Ln-dependent XoxF, source chain A.
- **Secondary only:** 6OC6, XoxF/La, source chain A. Its sequence is already
  represented by C5B120 in the calibration panel, so it is not an independent
  test.

The donor-bearing chain is unambiguous in all three structures. Each crystal
contains two enzyme copies. Chain A is frozen because the metal coordinates in
the already-used local 1H4I and 4MAE historical carves match the chain-A site
exactly. Chain A is also frozen for 6OC6 for consistency. Selection must never
fall back to “first ligand” or search all chains at execution time.

## Exact fixed cores

| Structure | Metal / PQQ | Protein fragments | Typed source CN at 3.1 A |
|---|---|---|---:|
| 1H4I | `A:CA701/CA`; `A:PQQ601` | `A:GLU177;A:ASN261;A:ASP303;A:ARG331` | 6 |
| 4MAE | `A:CE601/CE`; `A:PQQ602` | `A:GLU172;A:ASN256;A:ASP299;A:ASP301;A:ARG326` | 9 |
| 6OC6 | `A:LA701/LA`; `A:PQQ702` | `A:GLU192;A:ASN276;A:ASP318;A:ASP320;A:ARG345` | 9 |

The homologous +2 residue in 1H4I is `A:ALA305`; it is recorded but excluded
because it is not an acidic direct donor. The +2 residues in 4MAE and 6OC6 are
Asp and are included. The Arg is a required second-shell fragment in each case,
not a metal donor. Its `NH2` contact to catalytic-Asp `OD2` is 2.989, 3.082,
and 3.034 A in 1H4I, 4MAE, and 6OC6, respectively.

All PQQ residues have the canonical 24-heavy-atom PDB CCD naming graph and are
frozen as `pdb_ccd_pqq_v1`. The complete atom-level 3.1 A donor ledgers are in
`holdout_spec.tsv`; CN is descriptive and cannot change fragment membership.

## Waters, alternate conformers, and the 4MAE adduct

1H4I contains no source water or alternate conformer. 4MAE contains 1,085
waters, but none is within 3.6 A of chain-A Ce; the nearest is
`A:HOH796/O` at 5.091 A. Its four alternate-conformer residues are remote from
the selected core (`A:ASP162`, `B:ASP205`, `B:ASP207`, and `B:LYS317`). 6OC6
contains 26 waters, none within 3.6 A of chain-A La; the nearest is
`A:HOH807/O` at 9.133 A, and it has no alternate conformer.

4MAE has one material complication: `A:15P603/OXT`, part of a PEG-1500
crystallization adduct, directly coordinates Ce at 2.747 A. It is not part of
the conserved PQQ-MDH core. Removing it leaves a ligand vacancy in a fixed
crystal geometry. This cannot be hidden as routine cleanup: the primary 4MAE
result must be described as a **dry fixed-coordinate structural-transfer test**,
not as the energy of the intact crystallographic first shell.

## Mandatory dry policy at release

If—and only if—the 25-protein calibration gate passes, holdout preparation must
enforce the following policy and record it in the carve manifest:

1. Start from the exact hashed source in `holdout_spec.tsv` and select model 1,
   chain A only.
2. Retain the metal coordinate, complete PQQ, and exactly the listed protein
   fragments. Preserve all source heavy-atom coordinates.
3. Exclude every crystallographic water and every noncore heterogen. For 4MAE,
   this explicitly removes all of `A:15P603`, including its Ce-bound OXT atom.
4. Add no synthetic water, counterion, or replacement ligand; do not optimize,
   minimize, or fill the 4MAE vacancy.
5. Build paired La/Ca arms with byte-identical nonmetal coordinates. The arms
   may differ only in metal identity, charge, and electron count.
6. Fail closed on a source-hash, selector, residue-identity, altloc, atom-name,
   coordinate, charge-ledger, or paired-coordinate mismatch.

Any wet, adduct-retaining, vacancy-filled, or relaxed 4MAE calculation is a
different scientific experiment and requires a new protocol ID. It cannot
replace or rescue the preregistered primary result.

## Historical protonated files are not v3 inputs

Existing protonated 1H4I and 4MAE files are recorded in `holdout_spec.tsv` only
for provenance. Their ligand chains/residue numbers were remapped: the selected
1H4I metal/PQQ become `B:CA2` and `B:PQQ1`; the selected 4MAE metal/PQQ/15P
become `B:CE1`, `B:PQQ2`, and `B:15P3`. Historical preparation also reset the
4MAE Ce occupancy from 0.60 to 1.00. A future v3 preparation must therefore
rebuild from the hashed raw source under the pinned workflow and emit a fresh
manifest; it must not silently reuse either historical artifact.

## Interpretation boundary

These are structure-level transfer controls, not blind sequence validation.
The mature 1H4I sequence occurs within calibration protein P16027, prior
protocol results exist for 1H4I and 4MAE, and 6OC6 duplicates the C5B120
sequence family member. None may be used to revise the calibration chemistry,
threshold, or gap after energies are seen.
