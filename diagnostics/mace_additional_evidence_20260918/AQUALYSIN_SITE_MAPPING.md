# Aqualysin I: verified affinity direction, unresolved structural assignment

2026-09-18. Approved evidence curation only. No preparation, modelling, energies, jobs, new benchmark labels, or changes to earlier records.

## Decision

**Retain one direct, conditional Ca-preference observation. Do not promote either 4DZT calcium to a confirmed site-resolved control.** The strongest published hypothesis points to **4DZT author chain A, Ca303**, but does not experimentally identify it as the site measured in the La/Ca assay. The other site and the covalent inhibitor also require explicit treatment before any calculation.

## What the experiment establishes

Suefuji et al. report association constants **Ca²⁺ 6200 M⁻¹; La³⁺ 1100 M⁻¹**, using competition with sodium monitored by ²³Na NMR (Table 1, p.1285; model Eq.10, p.1283). Conditions: 22°C, pH6.0, 20mM MES-NaOH, 10% D₂O, 0.10mM holoenzyme; measured background sodium 6.53mM. No uncertainties accompany these constants. The linewidth model assumes bound competitor is negligible relative to free competitor. [Primary experiment, DOI10.1271/bbb.66.1281](https://doi.org/10.1271/bbb.66.1281).

The enzyme preparation retains one tightly bound Ca; the weaker Ca site is initially Ca-depleted. Thus this is a **solution binding-affinity direction at the weak site, conditional on the retained strong-site Ca**, not an activity or physiological-use label. The paper additionally identifies a Ca-insensitive La site, X; its approximate La association constant of 20000 M⁻¹ belongs to a different observation. Thermostabilization is not an affinity ranking.

The assay uses recombinant *Thermus aquaticus* YT-1 aqualysin expressed in *E. coli* from a precursor lacking its C-terminal prosequence, matured at 70°C. No point mutation or PMSF treatment is specified. Its methods reference a 1997 mutagenesis paper for preparation; that citation alone does not establish that the NMR sample was an Asn219 mutant. No sample sequence is supplied, so exact sequence equivalence cannot be independently pinned from the assay paper.

## Residue crosswalk: names alone are unsafe

These are **deposited connections**, not a new donor-selection calculation. Author atom/residue identifiers are used throughout.

| 4DZT ion | Deposited protein donors | Deposited linked water | Published nomenclature |
| --- | --- | --- | --- |
| A:Ca302; label chain C | Asp11 OD1/O; Asp14 OD1; Gln15 OE1; Ser21 OG; Ser23 O | A:HOH1012 O | Corresponds to VPR Ca3; called AQN Ca1 in the 2018 paper |
| A:Ca303; label chain D | Val170 O; Ala173 O; Thr175 OG1; Asp196 OD1 | A:HOH1074 O | Corresponds to VPR Ca1; called AQN Ca2 in the 2018 paper |

The [2005 VPR structure paper, calcium-binding-sites section](https://doi.org/10.1111/j.1742-4658.2005.04523.x) defines Ca1 around Pro171/Gly173/Asp196 and Ca3 in the N-terminal Asp/Gln loop. The [2018 primary comparison, Table7](https://pmc.ncbi.nlm.nih.gov/articles/PMC9085296/) lists AQN Ca2 as V174/A177/T179/D200 using its alignment numbering; these must not be mistaken for native 4DZT residue numbers. Its electrostatic simulations provide no independent experimental affinity label.

The [2019 VPR study, pp.160–161](https://doi.org/10.1016/j.bbapap.2018.11.010) argues that the aqualysin weak site likely corresponds to VPR Ca1, based on structural analogy and stability behavior. This traces to **4DZT A:Ca303**. The authors explicitly make an inference; they do not perform a site-resolved aqualysin La/Ca experiment. Their discussion also notes that homologous sites can have substantially different affinities. Full paper is preserved in the [author's dissertation](https://opinvisindi.is/bitstreams/60f398cb-139f-4797-851a-899a5d0245c9/download).

**Missing link:** a residue-specific assignment connecting the weak-site competition assay to Ca303, such as an existing selective-depletion structure or site-specific binding study. This bounded search did not locate one or a better matched deposited aqualysin structure. This is not proof that none exists.

## Structure/assay compatibility and preparation requirements

[4DZT](https://www.rcsb.org/structure/4DZT) is a 1.95Å natural-source YT-1 mature monomer: chainA residues1–276 map to P08594 residues128–403, with no reported sequence mutation. All276 residues are modelled. The current download is byte-identical to the previously archived source. Crystal conditions are pH7.5,296K,21% PEG800 and100mM phosphate. The 2002 assay uses different buffer/pH and recombinant production.

- **Covalent inhibitor:** PMS A301 is bonded to catalytic Ser222 OG (deposited `_struct_conn` and modification record). No such treatment is reported for the NMR sample. Retaining the adduct models an inhibited structure; removing it is a chemical preparation change requiring an agreed policy. It is not an innocuous solvent deletion.
- **Two calcium ions:** both deposited occupancies are1.00. A future conditional-site model must retain the other Ca and explicitly state that assigning it as the strong site depends on the same unconfirmed mapping. Swapping both metals changes the experimental question.
- **Waters and chemical state:** the linked waters above are source inventory, not an approved water-selection rule. Freeze protonation, disulfides, water selection, and the treatment of the inhibitor and spectator Ca before any calculation. Alternate conformers occur at Asp17 and His70; selection must be recorded. No preparation was done here.
- **Evidence use:** keep one biological group, `aqualysin_I_Taq_YT1`. Preserve the weak-site label and the unassigned La-only site separately. Neither two ion positions nor multiple structural hypotheses creates additional independent controls. Scores have not been inspected; the evidence/mapping has been consumed during development.

**Next decision:** use this as an affinity-evidence record with unresolved site mapping, or separately agree to an explicitly hypothesis-labelled structural test. It is not ready for a confirmatory site-classification denominator. Do not choose a site using its predicted score.

## Artifacts

- [Machine-readable result](AQUALYSIN_SITE_MAPPING.json): evidence, candidate selectors, limitations and SHA256 pins.
- [Source inventory](../../workspaces/mace_additional_evidence_20260918/aqualysin/4DZT_inventory.json): exact source atoms, deposited bonds, alternates and modification records.
- [Acquisition receipts](../../workspaces/mace_additional_evidence_20260918/aqualysin/acquisition.json) and pinned primary captures in the same workspace. The new J-STAGE download failed DNS resolution; the existing real six-page PDF was reused and pinned instead.
