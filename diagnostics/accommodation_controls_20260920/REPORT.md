# Matched accommodation controls — 2026-09-20

**A real Ca/La crystal pair is available: bovine alpha-lactalbumin 1F6S/6IP9. An explicitly matched two-water version already has all four DFT endpoints. Its La-source geometry shifts the Ca−La contrast by −10.4374 kcal/mol, toward Ca.** This is a conditional geometry/hydration result on consumed cases, not improved affinity discrimination.

No new electronic calculations, folds or jobs were launched. Baseline and existing jobs remain unchanged. The parent explicitly extended the initial curation task to this four-energy archive comparison after the equal-composition inputs were identified.

## Small useful set

| Control | Available inputs | Legitimate test | Limitation |
|---|---|---|---|
| Alpha-lactalbumin native pair | 1F6S chain A, Ca201; 6IP9 chain A, La203 | Observed accommodation within the same protein | Native cores have 2 versus 3 waters; different crystal packing/conditions |
| Alpha equal two-water intervention | Archived `1F6S__full` and `6IP9__minus_A_322`, 40 atoms each | Same-composition electronic contrast across source geometries | Deleted-water 6IP9 is an explicit intervention, not its native state |
| GGR source-geometry control | Existing 1GLG, 2FW0, 2FVY, 52-atom cores, no waters | Same-metal structural/context variability | One affinity group; sugar, packing and radiation differences confound causal attribution |
| XoxF1 D320A perturbation | WT 6OC6; experimental WT/mutant metal loading | Mutation-dependent loading as a future model test | No observed D320A or Ca-bound geometry in the paper; no mutant prepared here |

Exact source, preparation and endpoint paths/hashes are in [MANIFEST.json](MANIFEST.json). These are two ready structural control families and one mutation candidate, not four independent biological observations.

## Alpha atom mapping and chemistry

Both depositions specify the identical 123-residue P00711 sequence. Observed chain A contains 122 residues in 1F6S and 120 in 6IP9; unresolved tail atoms cannot be invented to make whole-protein compositions equal. The retained core is complete in both. The primary paper explicitly compares these structures and describes different packing and hydration. [La structure and analysis](https://doi.org/10.1038/s41598-018-38024-1), [Ca structure](https://doi.org/10.1074/jbc.M004752200).

The shared core contains 18 observed protein heavy atoms: Lys79 C/O with Phe80 N; Asp84 C/O with Leu85 N; and CB/CG/OD1/OD2 from Asp82, Asp87 and Asp88. Including the prepared protein hydrogens and seven caps gives 33 nonmetal atoms. These atoms map exactly by source identity; Ca/La paired coordinates are identical within each preparation. Caps and hydrogens are prepared coordinates, not crystallographic observations. Carboxylate OD1/OD2 naming does not imply unique physical oxygen identity.

An all-common-Cα alignment uses 120 atoms, with 0.43166 Å RMSD. The source water oxygen mapping is:

| 1F6S | 6IP9 | Separation after protein alignment |
|---|---|---:|
| A:HOH211/O | A:HOH326/O | 0.47741 Å |
| A:HOH212/O | A:HOH310/O | 0.27881 Å |
| No matched water | A:HOH322/O | — |

This is a positional correspondence, not the identity of persistent water molecules. It was determined without scores. The existing hydration-square experiment already prepared and scored **every** single-water deletion, including A322; it is not fresh validation. Both native inventories remain in the inventory. The paper describes eight-coordinate La, while the historical typed-distance check also counts Asp82 OD2 at 2.885 Å; no donor or threshold was changed here.

Mappings: `workspaces/accommodation_controls_20260920/alpha_atom_mapping.json`, `alpha_common_heavy.tsv`, `alpha_equal_two_water_pair.json`.

## Reused DFT matrix

All four exact XYZ hashes match successful execution receipts and normally terminated outputs. Inputs share ORCA 6.1.1 native r2SCAN-3c/CPCM(Water)/DefGrid3, archived NormalSCF default, Ca −1/La 0 singlets, and `reference_internal_geometry` water-H normalization. No contextual-water-prepared energy was substituted.

| Two-water geometry | E_Ca / Eh | E_La / Eh |
|---|---:|---:|
| 1F6S full | −1855.780459462680 | −1209.692433436387 |
| 6IP9 minus A322 | −1855.737201409227 | −1209.632542288746 |

With `R=E_Ca−E_La`, `ΔR=R(6IP9−A322)−R(1F6S)=−10.437424184904337112 kcal/mol`, converting once with 627.509474 kcal/mol/Eh. The two endpoint geometry changes are +27.14483837 kcal/mol for Ca and +37.58226255 for La. Both metals favor the 1F6S geometry within this restricted representation, with a larger penalty for La. There is no compatible affinity threshold or population estimate in this comparison.

This cannot establish that metal accommodation is absent or that the Ca geometry is universally preferable: hydration removal, fixed prepared H orientations, fragment boundaries, crystal environments and electronic treatment remain material factors. It does show that simply choosing the experimentally La-bound protein coordinates does not automatically improve the isolated-core contrast.

Full decimals, checks and receipts: `workspaces/accommodation_controls_20260920/alpha_archived_dft_matrix.json`. The four reused endpoints originally took 371.932646 summed endpoint wall-seconds on 16 ranks each (5950.922336 assigned rank-seconds); these are historical costs. New endpoint count and GPU usage are zero.

## Context construction options for the parent

The existing local contexts are **not** immediately compatible: they have metal-dependent water H coordinates, and source membership differs. 1F6S includes the peptide unit around Asp82–Asp83 absent from the 6IP9 expansion. After deleting A322, counts would be 106 versus 101, not equal composition.

The implementation-ready option is to freeze the union of chemically complete source peptide units, verify corresponding source atoms, rematerialize caps through the source graph, and copy every surviving core atom—including archived water H—from these exact endpoints. Remove A322 explicitly. Preserve the standalone core separately for any subtractive expression; a restored covalent context replaces synthetic caps. Do not use old context energies as though they belonged to the four-endpoint matrix. Exact source maps and required fixes: `workspaces/accommodation_controls_20260920/context_construction_options.json`. No new context or optimization was generated here.

## Mutation evidence and independent-group readiness

Good et al. tested XoxF1 WT/D320A metal loading in medium containing 20 µM Ca, with or without 2 µM La. With La, WT was 39% La-loaded with no detectable Ca; D320A instead contained Ca and only trace La. Without added La, WT was >97% Ca-loaded. This supports altered **loading**, not an equilibrium affinity constant; Ca loading did not restore mutant catalytic activity. The deposited crystals are WT La/PQQ 6OC6 and WT La without PQQ 6OC5. ExaF D319S supplies functional evidence only. [Primary paper](https://doi.org/10.1074/jbc.RA120.013227).

Existing 6OC6 preparation failed a zero-heavy-reconstruction rule: 130 missing heavy atoms across 44 noncore residues. The failed whole-protein preparation must not be reused as valid; an explicit local policy would be needed. No mutation coordinates were manufactured.

The existing ledger supplies **no additional fully ready, uniquely mapped direct La/Ca affinity group**. Aequorin has an archived ordered three-site vector but the competitive site is unmapped. Aqualysin has direct Ca>La evidence but unresolved structural site assignment. LanP has strong protein-level La/Ca competition evidence but coupled adjacent metals and an unresolved site-population mapping. PqqT has direct affinity data but unresolved metal placement. These limitations are retained in `independent_control_readiness.json` rather than promoted into labels.

A bounded extra lead, cod parvalbumin 9B26 versus 2MBX, is not a clean pair: La X-ray versus Ca NMR ensemble, G/M N-terminal difference, extra partially occupied La and crystal anions, and no verified matched La/Ca affinity. Carp parvalbumin labels cannot be borrowed. [La primary study](https://doi.org/10.1002/pro.5226), [Ca NMR primary study](https://doi.org/10.1002/prot.24664).

## Reproduce the archive checks

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/accommodation_controls_20260920/build_inventory.py
```

This performs real source/paired-coordinate checks, rigid mapping, hash/receipt validation and direct energy algebra; it launches no scientific executable. It preserves original inputs and outputs. The next useful experiment requires a declared common context, not another unqualified isolated-core rescore.
