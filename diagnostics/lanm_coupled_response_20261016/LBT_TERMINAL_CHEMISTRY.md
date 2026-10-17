# LBT3 terminal chemistry: primary evidence remains insufficient

Bounded review of the consumed Hatanaka 2020 paper, its complete supplement and
the 7CCO/7CCN deposition records. No reserved-family outcomes, new chemical
preparation, numerical science or jobs were used. Existing models are unchanged.

## Paper and supplement

The main Experimental Procedures call the material synthetic peptides, but do
not identify a peptide supplier, synthesis resin, terminal protecting/capping
groups, a synthesis protocol, or an experimental terminal-charge measurement.
No explicit C-terminal amidation or N-terminal acetylation is documented there.
The MD section assigns negative Asp/Glu side chains and reports holo charge −2;
that is a **simulation-state declaration**, not proof of experimental termini.
[Hatanaka et al., Scientific Reports 2020, Experimental Procedures](https://pmc.ncbi.nlm.nih.gov/articles/PMC7656248/).

Supplementary Fig. S1, page 2, gives `FIDTNNDGWIEGDELLA` without terminal-group
notation. Fig. S5, page 6, says the peptide contains five acidic residues and
“exhibits a charge of -4 at pH 6.” Neither its caption nor the other supplement
pages specifies which terminal chemistry or protonation assignment produces
that value. Both pages were visually checked against the preserved PDF.
[Original supplement, Figs. S1/S5](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-020-76527-y/MediaObjects/41598_2020_76527_MOESM1_ESM.pdf).

The two statements therefore remain different primary-source declarations.
The peptide −4 statement alone does not establish amidation, and the holo −2 MD
statement does not uniquely establish a terminal microstate. Neither warrants
silently correcting the existing preparation.

## Official deposition metadata

Both official mmCIF records declare the ordinary 17-residue sequence, PHE1 through
ALA17, with `nstd_linkage=no`, `nstd_monomer=no` and chain modification flag `N`.
Their chemical component lists contain ordinary amino acids and La/Lu only.
Neither declares an amidated or acetylated terminal group, a separate cap, or a
modified-residue record. These are **absent modification declarations**, not an
experimental synthesis certificate proving free termini.
[7CCO official mmCIF](https://files.rcsb.org/download/7CCO.cif),
[7CCN official mmCIF](https://files.rcsb.org/download/7CCN.cif).

The relevant categories are `_entity` / `_entity_poly` (lines 147–157),
`_entity_poly_seq` (168–184), `_chem_comp` (204–214) and `_struct_asym`
(chain modification flag, line 287) in both preserved files. Generic chemical
component atom templates are not proof that all terminal atoms were present in
the deposited coordinate model. No formal charge was inferred from absent OXT.

## Existing model and consequence

`LBT3_PLAN.md` and `LBT3_PREPARE.py` explicitly choose free N-terminal NH3+ and
C-terminal COO−, with D3/D7/D13/E11/E14 deprotonated. Preparation 12469 preserves
the deposited 252 peptide atoms and adds H2/H3/OXT, producing 255 peptide atoms
plus one trivalent ion, peptide −5/holo −2. Those are declared modeling choices.
Missing deposited OXT or terminal H atoms is not evidence that the experimental
peptide was amidated or otherwise capped.

**Conclusion:** no experimentally documented terminal modification has been
established by this bounded evidence review. Keep the current free-terminus
model conditional, preserve the SI discrepancy, and do not relabel it as an
experimentally verified chemical state. Resolving the experimental chemistry
would require an explicit synthesis/terminal specification beyond these records.

## Exact local evidence

- Previous consumed evidence/state notes: `MONOMER_TARGETS.md`, `LBT3_PLAN.md`
  and `LBT3_STATIC_REVIEW.md` in this diagnostic.
- Primary XML: `lanthanide_series_capability/lbt3/raw/PMC7656248_fulltext.xml`,
  SHA256 `3ed402373989c5c5ad3ebaf3e5000eb926fd4a459a6d86bcfe544f00ac71d81d`.
- Primary SI: same directory, `41598_2020_76527_MOESM1_ESM.pdf`, SHA256
  `a21276ce24a352c9c1d1761fe594ef356a3c6f714ec3e0e5ce069c807fa9c21b`.
- Paper reading extracts, page renders and pins:
  `workspaces/lanm_coupled_response_20261016/lbt_terminal_chemistry_evidence_v1/PAPER_EVIDENCE.json`
  (`56d8bcec9100be19f4d07619bec5a8e371dd525a38b065c7e0bc894ddfe71c67`).
- Current official deposition downloads: the same evidence directory's
  `deposition/7CCO.cif` (`473cc166df369b085b1fe030fa03f4fc31ef8677b28f2385e4d74f6462549807`)
  and `deposition/7CCN.cif` (`dfe1316b8319ed1205bca731115b4e9e3d194183f6b33b50f51ced5afb92f7fb`).
  An independent reviewer checked the complete entity/component declarations.
