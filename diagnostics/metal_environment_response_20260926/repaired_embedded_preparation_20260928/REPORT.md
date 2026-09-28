# Repaired full-region source contracts ready

All three reviewed C-bound-H repairs now have separately versioned full-region A/B inputs compatible with `metal_environment_lady_embedded_prepare.py`. Nine actual source/mapping/state checks pass per source. No electronic endpoint, optimizer, field calculation or Slurm submission occurred in this preparation step.

**Use v2**: `workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/{Hans8DQ2,Hans8FNR,Mex8FNS}/INPUTS.json`.

| Source | Input SHA256 | Transferred real H |
|---|---|---:|
| Hans8DQ2 |85770b1d42be2ebe92051d6596206ad25a09011ca5e9c5a40b5a3e0ed5f2cec0|630|
| Hans8FNR |8aa2b5a94c56a3f4052ce172d3f990272f2b8c6ba8c98996a1418c69da90fbe2|630|
| Mex8FNS |4b3b98983af0cb5fcdbf7b8f52e35063308201ce9df777a30484acf5146b7930|597|

Each contract pins original inputs and separate reviewed admission. Actual region sizes195/198/209 and all atom identities, charges, physical spin/state metadata, water inventories and spectator assignments are retained. These source contracts preserve physical Dy multiplicity6; the downstream frozen-f adapter separately records effective multiplicity1/ECP55. No effective state was silently written into physical preparation metadata.

A exactly uses admitted full-source coordinates. B is regenerated from that A by the unchanged +2° Asp85/Hans or Asp86/Mex rotation, including HB2/HB3. All unaffected physical atoms remain fixed; B covalent length residuals are at most7.8e−15Å. Caps and their physical Jacobians are regenerated from unchanged heavy anchors and remain coordinate-identical to the originals. Both target metals receive exactly the same geometry for each configuration.

Environment rows retain old charge values, recipients and source ordering, but every moved exterior C-bound H gets its repaired source coordinate. A/B field bytes match within each source; no original-field/new-core mixture occurs. The local charge ledger's dipole diagnostics are recomputed at repaired coordinates, with original diagnostics separately retained. Charge closure remains the source's original QM/MM split. No missing ion/cross-LJ/full-hybrid terms are manufactured.

`physical_source_A.json` and `physical_source_B.json` retain exact full source IDs and coordinates; `core_mapping_A/B.json`, `environment_atoms.json`, `boundary_mapping.json` and the four La/Dy XYZs map directly to those physical arrays. Source protonation/parameter provenance remains original and is explicitly distinct from the new coordinate preparation. The original ff19SB export is not misrepresented as having repaired coordinates.

An initial v1 preparation produced identical new physical coordinates, endpoint XYZs and fields but retained original boundary dipole diagnostics. It is preserved unused; v2 fixes that diagnostic metadata, not chemistry. No existing input or active root electronic job was modified.

Reproduce into an unused directory:

```bash
python diagnostics/metal_environment_response_20260926/repaired_embedded_preparation_20260928/PREPARE_v2.py workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_reproduction
```

This runs the nine checks on each actual source and writes per-source `CHECK.json` plus aggregate `CHECKS.json`. It performs no molecular evaluation. Root owns any subsequent versioned electronic matrices; old energies cannot fill repaired cells. Improved H geometry is established, improved La/Dy discrimination is not yet established.
