# Alpha-lactalbumin: read-only failure diagnosis

Requested question: why does the discriminator miss alpha-lactalbumin? This audit inspected existing inputs, source structures, manifests, results and primary literature. It launched no scientific calculation and changed no preparation or label.

## Working diagnosis

The demonstrated problem is poor transfer of a frozen local-core electronic descriptor to a condition-qualified solution-affinity ordering. The leading physical hypotheses are incomplete hydration/reorganization treatment and loss of the surrounding protein's interactions in the carve. Neither has been isolated as the cause. Successful PQQ functional-class separation does not establish that the same descriptor orders relative affinity across different binding-site architectures.

## What the existing artifacts establish

- All ten new endpoints converged and passed receipt/algebra verification. Alpha ranks 14.559245/17.542705 kcal/mol below the previously computed Ca-favoring GGR using the same repaired protocol. The common reference cancels; a reference shift cannot fix that ordering.
- Both Ca-derived 1F6S and La-derived 6IP9 geometries were evaluated. Within either endpoint pair, coordinates, donor membership, protonation and water inventory are identical. Neither pair describes the free-energy cost of switching between different coordination/hydration states.
- Both alpha cores retain Asp82/Asp87/Asp88 side chains, Lys79/Asp84 carbonyls, and their respective two/three waters. No omitted directly coordinating residue was identified against the inspected structural description.
- The repaired graphs retain Lys79 C–Phe80 N and Asp84 C–Leu85 N. The original formaldehyde-like backbone defect is absent. Chemically valid caps do not establish that a fragment reproduces the intact protein's electrostatics.
- Alpha and GGR both have selected ligand charge −3, hence La-core charge 0 and Ca-core charge −1. Gross core charge alone does not explain their reversed ordering. Their donor types, waters and geometries differ.
- The alpha carbonyl fragments do not retain the entire parent residues: for example, the Asp84 side-chain charge is outside the selected core. This is a recorded limitation of the local representation, not an undisclosed new preparation defect or a demonstrated cause.
- 6IP9 La occupancy is 0.70; 1F6S Ca occupancy is 1.00. Selected donors/waters have occupancy 1.00 and no alternative conformers. Partial metal occupancy is a structural caveat, not a diagnosis.
- Our distance-based 6IP9 descriptor counts nine donors, including Asp82 OD2 at 2.885 Å. The paper calls the site eight-coordinate. The complete carboxylate is already retained either way, so this descriptor discrepancy does not change the calculated atom inventory.

## Primary evidence and its limits

The 1996 competition study supports La-over-Ca at the bovine strong site; its exact conditions remain unrecovered. That direction is retained. [Primary abstract](https://pubmed.ncbi.nlm.nih.gov/8652630/).

The 2019 study reports changed coordination/waters and La-only ITC with binding enthalpy +12.5±1.0 kcal/mol and entropy +64.4 cal/mol/K. It therefore provides a concrete example of favorable association supported by entropy. It supplies no measured **La-minus-Ca** entropy correction and cannot identify the missing term in our score. Crystal conditions (pH 6, 2 M sulfate), ITC (pH 7.4) and our preparation (pH 7) also differ. [Primary study, Figures 1–4 and Methods](https://doi.org/10.1038/s41598-018-38024-1).

Static electronic energies with continuum solvent are not experimental binding enthalpies. Do not equate the reported +12.5 value with an endpoint energy or add the reported entropy as a correction. Crystal B-factors do not provide a solution entropy estimate.

## Consequence

Retain alpha as a challenge case with its conditions qualification. We have not demonstrated a correctable coding defect, established a unique missing-physics explanation, or validated a replacement score. A next numerical diagnostic would require a specifically defined, matched preparation/energy cycle; changing waters, protonation, caps or geometry until the sign agrees is not a diagnosis. Existing failed global/environmental implementations remain archived.

Artifacts: `preparation_configs/{1F6S,6IP9}.json`, the two repaired manifests under `workspaces/benchmark_set_20260915/prepared/alacta_*_v2/`, `SCORING_RESULT_1199508.json`, and the pinned primary captures under `workspaces/benchmark_set_20260915/evidence/nonpqq/`.
