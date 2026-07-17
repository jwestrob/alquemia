# DFT Ca²⁺/Ln³⁺ Discriminator — Methods Paragraph for Colin's Paper

**Status:** Final wording (2026-07-13). Citations moved out of the paragraph by Colin; they will appear in the References section as noted below.

**Companion paper:** methodology paper in preparation (West-Roberts et al., in prep.) which will report full computational details, benchmarking, and validation. The paragraph below is deliberately hedged to reserve those details for the second paper.

---

## Final paragraph text (as of 2026-07-13)

For each protein co-folded with La³⁺, we constructed a quantum-mechanical cluster model of the metal-binding site by extracting the metal ion, first-shell coordinating residues, and inner-sphere waters from the predicted structure. Single-point density-functional calculations were performed on Ca²⁺- and La³⁺-loaded versions of each cluster in implicit aqueous solvation, holding the protein geometry fixed across the metal swap. ΔΔE (in kcal/mol) was defined as the cluster's Ca²⁺→La³⁺ swap energy minus the corresponding swap energy for the aquo ions in bulk water, isolating the site-specific contribution to metal preference. Positive ΔΔE indicates a site favoring Ln³⁺ over Ca²⁺ binding relative to bulk water. Full computational details, benchmarking, and validation will be reported separately.

---

## What this paragraph deliberately does NOT say

Preserved for Paper 2 (methodology paper):

- The composite functional (r²SCAN-3c)
- The specific DFT software (ORCA 6.1.1)
- The implicit solvation model (CPCM)
- The basis set + ECP scheme (def2-TZVP + def2-ECP for La)
- The carve cutoff (3.0 Å for `carve_generic.py`, 3.2 Å for `carve_with_pqq.py`)
- The residue-type whitelist (Asp/Glu/Asn/Gln/Ser/Tyr — Cys/His/Met/Lys/Arg deliberately excluded)
- Backbone-O carving policy
- Aquo reference numerical value (DIFF_AQUO = −646.064555551745 Ha)
- Link-H capping and charge-assignment rules
- The asymmetric-carve outlier handling
- The CN/bidentate filters
- The benchmark validation against the Ca/Ln standards panel
- Cluster-size sensitivity results

## Terminology that MUST NOT appear

- **"Alchemical free energy perturbation"** — technically wrong AND misleading. We do not run FEP: no λ-windows, no MD sampling, no thermodynamic cycle. Our ΔΔE is a static electronic energy difference from single-point DFT. Any occurrence must be replaced.
- **"Free energy"** — ΔΔE is an *electronic* energy, not a free energy (no ZPE, no entropy, no finite-T).
- **"12-6-4 model"**, **"polarizable force field"**, **"MD sampling"**, **"potential of mean force"** — not applicable, we don't do MD.

## Recommended replacements for prior FEP-framed text

Any place in the abstract or body that previously said "alchemical free energy perturbation" should become **"quantum-chemical metal-selectivity scoring"** or **"cluster-based DFT metal-selectivity scoring"**.

Any place that said "predicted free-energy difference" should become **"predicted energy difference"** (drop the word "free").

## Citations (moved to References section by Colin)

See `citations.bib` in this directory for BibTeX. Six references cover the software stack this paragraph implies:

| # | Ref | Purpose |
|---|-----|---------|
| 1 | Neese 2022 WIREs Comput. Mol. Sci. 12, e1606 | ORCA 6.1.1 |
| 2 | Grimme et al. 2021 J. Chem. Phys. 154, 064103 | r²SCAN-3c composite method |
| 3 | Weigend & Ahlrichs 2005 PCCP 7, 3297 | def2 basis + ECP nomenclature |
| 4 | Barone & Cossi 1998 J. Phys. Chem. A 102, 1995 | CPCM implicit solvation |
| 5 | Cao & Dolg 2004 J. Mol. Struct. THEOCHEM 673, 203 | Small-core Ln ECP |
| 6 | Eastman et al. 2013 J. Chem. Theory Comput. 9, 461 | PDBFixer / OpenMM (protonation) |
