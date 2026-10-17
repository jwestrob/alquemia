# 6MI5 model 1 versus the current 8FNS model

## Decision

**6MI5 supplies a second Mex protein conformation, but it is not a matched second source for the current 2049-atom, four-metal, 167-water model.** The common protein sequence and EF3 heavy-atom identities map exactly. Construct, metal occupancy/species, and explicit hydration do not. A comparison retaining both deposited inventories would combine these state changes with conformation; it could not isolate a source-geometry effect.

This is a read-only inventory of existing deposits and preparation metadata. No coordinates were transformed, distances evaluated, preparations generated, or molecular calls submitted. No reserved outcomes were accessed.

## Deposited inventories

| Item | 6MI5 model 1 | 8FNS model 1/current source |
|---|---|---|
| Method / chain | NMR ensemble, selected model 1 of 12; author chain X | X-ray structure; author chain A |
| Construct | 117 represented residues, author 23–139 | 105 represented residues, author 29–133 |
| Sequence relationship | Native `PTTTTK`23–28, the exact common105 sequence below, then His6 at134–139 | Exactly the common105 sequence; no substitutions |
| Reference annotation | C5B164 residues23–133;134–139 explicitly annotated expression tag | C5B164 residues29–133 |
| Deposited ions | Three YT3: `X/200//Y`, `X/201//Y`, `X/202//Y`; occupancy1 each | Four Nd: `A/201//ND` through `A/204//ND`; A202 occupancy0.99, others1 |
| Deposited waters | None | 167 HOH oxygen records, author A301–A467 |
| Deposited hydrogens | 854 protein H records in model1 | None; current prepared model adds/provides H under its recorded preparation |

The common105 sequence (author29–133) is:

```text
VDIAAFDPDKDGTIDLKEALAAGSAAFDKLDPDKDGTLDAKELKGRVSEADLKKLDPDNDGTLDKKEYLAAVEAQFKAANPDNDGTIDARELASPAGSALVNLIR
```

The six N-terminal residues in6MI5 are native reference sequence, not part of the expression tag. Lack of deposited NMR waters is not evidence that the experimental protein has dry sites.

## EF3 and graph correspondence

The6MI5 deposited `struct_conn` records assign Y200 to the35–46 donor region (EF1), Y201 to59–70 (EF2), and Y202 to Asp84, Asp86, Asp88, Asp92 and Glu95 (EF3). These are deposited connectivity annotations; no new contact calculation was performed. There is no fourth deposited ion.

The accepted8FNS preparation identifies target EF3 as `A/203//ND`, with fixed Nd spectators at EF1/A201, EF2/A202 and EF4/A204. Its integer modeled occupancy includes A202 despite the deposited0.99 occupancy. It retains2045 protein/water atoms plus four ions, with protein/water charge−10 and whole formal charge+2. The local209-atom partition has204 real source atoms, one target and four artificial caps; its declared charge is0, with physical La(III) multiplicity1 / Dy(III) multiplicity6.

All **99 heavy protein atoms** selected in that209-atom core exist in6MI5 model1 with the same author residue number, atom name and element after the explicit chain X→A identity alias. The selected protein region includes partial residues66/67 and83/96, and residues84–95. Cap retained/omitted anchors are C66/CA66, N67/CA67, C83/CA83 and N96/CA96. The Asp84/Asp86/Glu95 CA–CB and Leu83 N–CA coordinate anchors are therefore present at the identity level. This does not qualify new cap coordinates/Jacobians, hydrogen naming/protonation, contacts, or the full moving/fixed graph.

The exact current local region also contains seven whole waters: A320,331,358,373,376,394,434. None has a deposited6MI5 counterpart. The existing covalent export graph and physical/charge maps belong to8FNS; they are reusable as a declared topology specification only after a separate source-specific preparation and validation, not as already admitted6MI5 artifacts. Retaining6MI5's extra residues also changes full-graph membership, including a downstream hinge component.

## What a matched modeled state would require

1. **One common construct and terminal chemistry.** Matching the current105-residue model requires explicitly removing6MI5 residues23–28 and134–139 and preparing new termini at29/133. Those residues are internal in the deposited117-residue construct. Simply deleting records does not establish the existing protonation, hydrogen map or−10 charge. Keeping the full117 construct instead requires extending8FNS with chemistry absent from that deposit and re-preparing both sources.
2. **One declared metal-occupancy state.** To retain the current conditional EF3 swap with three Nd spectators,6MI5's EF1/EF2 Y ions must become modeled Nd spectators, EF3 Y must become the La/Dy target, and an EF4 spectator absent from6MI5 must be introduced with an explicit placement hypothesis. Three-site Y and four-site Nd states are not matched occupancy.
3. **One explicit hydration inventory.** The current167 waters, including all seven local waters, cannot be inherited from6MI5 observations. Any transfer or addition would be a declared modeled hydration construction requiring its own physical checks; a chain alias supplies no water correspondence. No such construction is proposed or performed here.
4. **Source-specific physical and boundary admission.** Rebuild/validate full atom IDs, terminal/protonation charges, covalent and protected/moving membership, caps and their Jacobians, field-charge shifts, FF export, and the shared energy-rule ledger before any comparison. Identical core heavy-atom names alone do not authorize reuse of the current209/2049 state or its force qualification.

A reduced common three-site or water-free model would likewise change the accepted8FNS state and Hamiltonian preparation; it would be a new declared comparison, not a second source under the current model. Thus6MI5 is feasible as a future **explicitly modeled source-state sensitivity study**, but the existing files do not support a deposit-preserving matched2049-atom comparison or a clean matched-source common pool.

## Exact existing evidence

- [6MI5 deposit](/groups/banfield/projects/environmental/sr/srvp2020/protenix/mmcif/6mi5.cif): SHA256 `9d646fad70cd65e5d4841cecae41a1ec07c61f4fd0ea466d46221e1e4a3cc31b`; `entity_poly`, `struct_ref_seq`, `struct_ref_seq_dif`, model1 `atom_site`, `struct_conn`.
- [8FNS deposit](../../workspaces/benchmark_set_20260915/evidence/lanthanide/8FNS.cif): `bb593dc407de9f1d00138ea19c417a503b88bae3a2926967d53595da312b2d9f`; sequence/reference and model1 atom records.
- [Accepted repaired INPUTS](../../workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Mex8FNS/INPUTS.json): `4b3b98983af0cb5fcdbf7b8f52e35063308201ce9df777a30484acf5146b7930`; charge, occupancy, water and source provenance.
- [Core mapping](../../workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Mex8FNS/core_mapping_A.json): `684e77b98d11faf7050b33fa07ed86f449b565ef0d81fc8f8a6c6f0bf435301d`; exact209 membership and cap anchors.
- [Boundary mapping](../../workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Mex8FNS/boundary_mapping.json): `a8cf4fff80bd0ad1472775ed30a3ef97bf9b307e2c844a1e836b9805be3d34e1`.
- [Existing covalent graph](../../workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Mex8FNS/export_v1/bonds.json): `7b1293cc72b7279aad323cb43d74d4014b769f1d986ccde46eb4797573395677`.
