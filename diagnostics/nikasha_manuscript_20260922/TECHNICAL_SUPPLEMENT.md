# Released fast PQQ protocol: technical supplement material

This describes the **already released static protocol**, not the September22
accommodation variants. The source of truth is `params/pqq_fast_v1.json`, SHA256
`5ccc66bba69ceddd59b9178f04618f5b65149f275f0a4b967c0213d29c71bfd7`.
Exact inputs and receipts below are from the completed 1H4I release operation;
the other27 reference preparations have their own recorded compositions/charges.

## Model and source preparation

- Public workflow identity: `fast_PQQ_OMOL_GFN2_ALPB_v1`.
- Energy identity: `native_OMOL_plus_native_GFN2_ALPB_transfer_v1`.
- Source preparation: `fixed_core_PQQ_source_to_complete_context_v1`; context
  `native_OMOL_complete_polar_context_disulfide_v2`. Existing canonical source
  normalizers, seeded standard-protein protonation, fragment policies and complete
  polar-context expansion are reused. Explicit source assembly, site/PQQ and
  five homologous residue-role selectors are required. The source request and
  full source-to-context mapping are preserved; a residue name is not a label.
- MACE checkpoint: `MACE-omol-0-extra-large-1024.model`, SHA256
  `9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`.
  Native `omol` head, float64, nonperiodic,6Å model cutoff, explicit total charge
  and multiplicity. No extra dispersion or solvent energy is added to MACE itself.
- Pinned software: mace-torch0.3.16, torch2.8.0, graph-longrange0.4.4. Complete
  lock and source inventory are linked by
  `workspaces/mace_omol_20260917/software_v1/software_manifest.json`, SHA256
  `7d2448b279b5191a529ef329a60442635601576794f9737efc64343c9a6e0168`.
- Ca and La have identical paired coordinates and nonmetal membership, with
  endpoint-specific element and total charge; both have multiplicity1. These
  dry PQQ preparations contain no explicit waters. No water reservoir, redox
  change, protonation ensemble, geometry optimization or entropy term is included.

## Solvent endpoints and actual input

ORCA6.1.1 native GFN2, executable SHA256
`38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83`.
Use matched vacuum/water-ALPB calculations at each metal's exact MACE geometry.
The actual released 1H4I Ca ALPB input is:

```text
! Native-GFN2-xTB NoAutostart ALPB(Water)
%maxcore 2000
%method
 WriteXTBParam true
 ReadXTBParam false
end
%scf
 SmearTemp 300
 UseXTBMixer true
end
* xyzfile -2 1 core.xyz
```

The exact vacuum input omits `ALPB(Water)`; the actual 1H4I La input uses
`xyzfile -1 1 core.xyz`. These charges are an example, not universal site charges.
Other states come from their prepared manifests. The runner inserts the recorded
MPI allocation and records the rendered input. `SmearTemp` is an electronic SCF
setting; this is not an equilibrium simulation at a declared biological temperature.
Native parameter exports, charge sanity, SCF convergence and normal termination
are checked. This released input has no research `MaxIter 500` override.

Actual solvent inputs/coordinates/output receipts reside under
`workspaces/pqq_fast_release_20260920/standard_1H4I_v1/prepared_score/scoring/cases/1H4I/solvent/`.
Its manifest pins all four endpoint inputs. Ca-coordinate SHA256 is
`d2f0908cd217272b11342369d14e488f2e0c351586cead815a4110963a0e4c18`;
La is`3056a1bf62c4f8f1e274043e68e66c22230171cb11ea54b48b9ab656e7f7009d`.
Different hashes include the intended element substitution; physical coordinates
are paired. Native MACE requests and results are pinned by the adjacent
`prepared_score/scoring/mace_manifest.json` and final standard result.

## Algebra and frozen decision policy

For each metal, `E = E_OMOL,vac + E_GFN2,ALPB − E_GFN2,vac`.
Then `R = E_Ca − E_La`. Convert MACE eV by23.06054783061903 and GFN2 Hartree by
627.509474, each exactly once. The latter is the retained released calibration
factor. Store unrounded components before reporting a rounded contrast.

The published context bands are:

| Decision | R in this protocol's model kcal/mol |
|---|---:|
| Ca-supported | ≤−405464.18774828 |
| La-supported | ≥−405459.1113819997 |
| Inconclusive | strictly between the two bands |
| Unavailable | missing/unsupported/nonconverged required component |

These bands derive from the25 designated canonical development inputs in
`workspaces/compact_solvation_20260920/full_v1/comparison_v1.json`, SHA256
`0fbbf688f922ae8e1c32d26df79284de6897c406cd11465dfed302fd48dd48a0`.
There is no compatible absolute aquo reference for this composite, no universal
zero threshold and no affinity/probability interpretation. Never import its bands
into a changed geometry/state protocol without an explicit transfer test.

## Actual release and comparison artifacts

`diagnostics/pqq_fast_release_20260920/RELEASE_PANEL.json` pins the28 reference
preparations and measurements. The complete standard CLI output is
`workspaces/pqq_fast_release_20260920/standard_1H4I_v1/result_1203729.json`.
Its `plan.json` pins executable implementations and source request. The measured
release used oneH200/32CPU/200000MiB per standard allocation, about42.3seconds
per source excluding folding. Two MACE and four GFN2 calls are needed per source.
The additional structural comparison is described in the
[methods outline](METHODS_OUTLINE.md); repeated folds are not independent labels.

All later shared-pool/angular/joint results and numerical iteration qualifications
have separate identities. None changes this supplement's released recipe.
