# Nine whole-chain occupancy states are ready

All three declared sources support their complete chain A, all chain-associated
crystal waters and all three occupancy hypotheses. The preparation emits **nine
states and 18 paired La/Dy XYZ files**, without missing heavy-atom reconstruction,
peptide caps, new protonation, optimization, energies or Slurm jobs. This enables
the parent's whole-scaffold pilot; it does not establish within-series accuracy.

| Source | Protein residues/atoms | Crystal waters | EF12 or EF23 atoms/charge | EF1234 atoms/charge |
|---|---:|---:|---:|---:|
| Hans 8DQ2 A | 110 / 1665 | 74 | 1889 / +2 | 1891 / +8 |
| Hans 8FNR A | 110 / 1665 | 179 | 2204 / +2 | 2206 / +8 |
| Mex 8FNS A | 105 / 1544 | 167 | 2047 / −4 | 2049 / +2 |

The integral protein charges are −4 for Hans and −10 for Mex under the existing
pH5 protonation. Waters are neutral; each selected ion adds +3. La and Dy share
identical coordinates, atom order, H inventory and total charge. No charge or
protonation adjustment follows occupancy or endpoint metal.

## Exact delivered interface

`workspaces/lanm_global_occupancy_20260923/prepared_v1/manifest.json`

SHA256: `2b0b906300fe24b1147813885c32c02253d26d1f50b974e316a52fc1fcf342f5`.

- `sources/<source_id>/source.json`: experimental/protonated-source receipts,
  sequence length, residue templates, charge, complete-chain connectivity,
  site identities and water inventory.
- `sources/<source_id>/atoms.json`: all protein/water atoms and covalent bonds;
  `audit.json` retains the original H coordinates and every radial H adjustment.
- `states/<source_id>__<EF12|EF23|EF1234>/state.json`: atom count, metal indices,
  occupied/removed sites, common mapping, La/Dy XYZ/state pins and source receipt.
- Each `mapping.json` has `atoms[]` with zero-based index, physical/source IDs,
  chain, residue/name, element, source and prepared coordinates. `bonds[]` holds
  actual protein/water covalent index pairs, identities and lengths. **No metal
  coordination bond or force constant is imposed.** Protein and waters precede
  the explicitly indexed ions; downstream code must use `metal_indices`.

## Source completeness and preparation

The original heavy-atom inventory agrees with the retained pH5 preparations.
Complete peptide connectivity is supported: 109 peptide bonds for each Hans
chain and 104 for Mex, no missing internal atoms or terminal completion. The
8FNR protonation receipt's missing-residue annotation concerns another chain;
selected chain A is complete. All standard residue templates match ff19SB.

The reused protonation files contain long generated X–H bonds (original protein
ranges 1.155–1.305 Å across sources). The already established ff19SB radial-H
projection restores each templated bond length, preserving attachment, direction
and number. Retained water H directions use O–H 0.9572 Å; their inherited H–O–H
angles span 103.30–106.90°, without angle relaxation. Maximum H displacements are
0.2404, 0.2699 and 0.2392 Å for 8DQ2, 8FNR and 8FNS respectively. This is declared
input preparation, not an optimized proton ensemble. Every H has one covalent
parent. No all-atom pair is closer than 0.5 Å after preparation.

OpenMM Å↔nm serialization introduces at most **1.63×10⁻¹⁴ Å** heavy-coordinate
roundoff; the explicit numerical preservation check is 10⁻¹² Å, with zero physical
heavy movement. An initial exact-floating-list assertion failed on this roundoff;
the test alone was corrected, with `TESTS_v1.txt` retained. Prepared files and
implementation were not changed or rerun. The final five real-fixture checks pass
in 3.304 s and include full source/graph, fixed water/H inventories, paired XYZ,
charge/spin/parameter parity and actual EF4 identity. Preparation took 14.090 s
locally; no allocated compute or molecular calls.

## Occupancy and state interpretation

EF1–4 map to deposited author residues A201–204. 8DQ2 has La201–203 and Na204;
its four-Ln state explicitly replaces Na at EF4. 8FNR has deposited Dy201–204,
and 8FNS has Nd201–204. Two-ion states remove the nonselected deposited ions.
These hypotheses do not assert equilibrium occupancy, nor infer folding from
four-site loading. Each source retains its own full author-chain water inventory
across all states; counts intentionally differ between crystals.

Physical La multiplicity is 1. Dy uses the predeclared maximum-spin hypotheses
11 for two ions and 21 for four; other inter-ion couplings and spin–orbit states
are untested. Actual pinned native GFN2 parameters have 5d/6s/6p shells with
three reference valence electrons for both metals. Its distinct f-in-core
representation uses an effective singlet. All-electron/physical spin parity and
native valence parity pass; native counts are 5264 (8DQ2), 6104 (8FNR) and 5690
(8FNS), equal across occupancies because extra Ln(III) contributes no net valence
electrons. This does not validate a physical magnetic ground state.

Only the complete **monomer chain A** is represented. Hans metal-coupled dimer
and folding equilibria, missing crystal environment, protonation changes and
water exchange are outside this state model. Same-occupancy Hans/Mex exchange
comparisons can cancel elemental/aqueous offsets; two-versus-four total energies
are not occupancy free energies. All source/occupancy results must remain visible.
