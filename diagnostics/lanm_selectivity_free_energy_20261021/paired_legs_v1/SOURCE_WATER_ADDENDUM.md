# Deposited inner-water inventory

The pinned 9CEQ source is a **solution-NMR ensemble**, with 20 models and 5,160 atom records. Every model contains 257 peptide atoms and one La ion. There are no deposited water atoms or water entities in any model. In model 1, the sole nonpeptide atom is LA in chainB, atom 258, occupancy 1.00; all 258 model 1 atoms have occupancy 1.00. No source metal–water distance or water occupancy can be reported because no water is modeled.

The original `PREPARE_PROTEIN_V2.py` asserts exactly 258 source model 1 atoms and 257 peptide atoms, then retains the peptide and metal separately. **It did not discard deposited coordinating waters.** The complete model-by-model labels, components, occupancies and empty water-distance inventories are preserved in `workspaces/lanm_selectivity_free_energy_20261021/paired_legs_v1/SOURCE_WATER_AUDIT_v1.json` (SHA256 `a6c1e9f3cb3cdb5c979b55ad96e690d2670c45452a598f1a1235ae279cf7a1e0`).

Absence of waters in a modeled NMR ensemble is not evidence for a physical hydration number q=0. It also does not establish q=1. **Inner-sphere hydration is unspecified by these coordinates.** The bound packing's closest water-O distance of5.06324 Å is therefore a preparation state, not an experimentally supported hydration state.

Before a production exchange calculation, hydration-state preparation and adequate exchange sampling must explicitly resolve this missing population. A short stability pilot or hoped-for water ingress cannot establish equilibrium occupancy. No hydration number was selected from classifier outcomes; no geometry, parameter, or native calculation was changed for this inventory.
