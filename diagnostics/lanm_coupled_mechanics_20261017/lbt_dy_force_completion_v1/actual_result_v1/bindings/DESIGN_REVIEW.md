# Independent design review

**Scientifically coherent and informative.** Completing the missing force arrays tests whether the reported Dy energy agreement also describes local forces, and whether model/reference force errors cancel in a matched Dy-minus-La contrast. Three fixed calls are a contained way to fill this evidence gap. The comparison does not require the separate step725 La reference.

Reviewed SPEC SHA: `5438ca99d6166e7164fe980e5461614254d6c4e2185712ddcf6fa08db31071ed`.

## Actual evidence and denominator

I read the complete collectors `9229f7be` (7CCO) and `fff95c2f` (7CCN), their reports, and all eight pinned full256 force tables. The tables confirm five available MACE comparisons and exactly these three missing arrays: 7CCO/before_Dy, 7CCN/before_Dy and 7CCN/after_Dy. All eight native reference gradients are available.

Final inventory `e2c196d1b5961cc74d50647e0565780928c3c1c103b282d563833f4b7e003957` preserves this denominator and each false force/state qualification flag. Its bounded cache check covers 16 unique POINT pins in the two already-bound complete pool menus; exact matching missing-cell entries are also energy-only. That supports the three-request scope within the declared source inventory, without claiming an exhaustive filesystem search. The inventory's model-checkpoint identity is inherited from the qualified configuration; runtime verification remains required.

Retain both original before/after pairs: the 7CCO hybrid-La deformation and the 7CCN whole-MACE-La step80 deformation. These are two consumed deformations of the same peptide. Do not substitute the later model-converged geometry or interpret them as independently equilibrated La/Dy endpoints.

## Required invariants

- Freeze all eight cells, three-or-fewer missing requests, five reused force authorities, per-metal solvent receipts and native reference provenance before model loading. Reuse any exact compatible force supplement found by inventory.
- Require identical coordinates and physical-ID order within each matched La/Dy geometry, with only the declared target/state substitution. Build `F_Dy − F_La` separately for each model, then compare those vector differences. Keep the four matched comparisons and eight endpoint slots explicit when any component is unavailable.
- Preserve physical Dy multiplicity6 for MACE, charge−2, whole256 atoms, float64, original model/head/helper/runtime and force sign. Never substitute the DFT effective multiplicity1 into MACE to make labels appear matched.
- The existing `1e−7 eV` energy-replay tolerance is present in original whole-force configuration `d3c28300`. It is a reproducibility check for each new force request, not a force-error or physical-accuracy threshold. Preserve original energy-only records and reported energy/work comparisons; append supplements separately.
- Admit finite full forces, all raw input/state/provenance checks and saved comparisons before/around the bounded calls. Stop scheduling on the first admission failure; preserve partial/unavailable rows without retries or tolerance adjustment.

The two collector schemas differ. In particular, 7CCO force tables expose whole-MACE composite forces but do not expose a whole-MACE vacuum column, while 7CCN tables expose both. Recover reused vacuum forces from their pinned native records. Preserve each source's established gradient-conversion arithmetic; do not introduce a cross-order bitwise replay requirement or silently loosen a gate to accommodate roundoff.

## Interpretation boundary for Dy

Physical-multiplicity6 MACE and spin-free frozen-4f9/ECP55 effective-singlet DFT are **different approximate electronic representations** of the declared DyIII system. Their comparison is useful for this conditional mechanics question, but it is not a spin-matched ground-state validation. A Dy-specific discrepancy cannot uniquely identify a MACE defect, frozen-core error or omitted spin-orbit contribution. Good metal-error cancellation would likewise not validate physical selectivity. The SPEC's inherited state caveats should remain attached to the new matched-metal results, not only the individual endpoint records.

Shared ALPB cancels from each model-versus-reference force error and its Dy-minus-La difference. It need not cancel from each model's Dy-minus-La force vector itself because the two metals can have different solvent terms. Preserve each metal's own cached solvent; no independent solvent validation follows.

Full vector norms/errors and original frozen projections answer the declared question descriptively. No new agreement threshold is needed. Nonstationary endpoints, conditional termini/protonation, vacuum-dianion assumptions and nonuniform native SCF/gradient limitations remain. No affinity, physical minimum, optimizer validation or natural-LanM transfer follows from completing this panel.

Review used source and saved metadata only. Inventory/config/source readiness and the actual bounded execution remain separate gates; no scientific imports, numerical replay, jobs or production edits occurred.
