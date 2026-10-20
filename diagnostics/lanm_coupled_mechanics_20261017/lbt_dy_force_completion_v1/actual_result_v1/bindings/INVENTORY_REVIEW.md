# Three missing Dy force requests: input review

**Confirmed: eight reference cells, five existing MACE force arrays, three energy-only MACE requests.** Completing the three fixed requests is independent of running La reference15545. It fills a meaningful force-comparison gap; coverage alone cannot establish physical force accuracy.

## Exact denominator

| Source | Geometry | La MACE force | Dy MACE force | Whole PBE0 gradients | Fixed-FF ALPB gradients |
| --- | --- | --- | --- | --- | --- |
| 7CCO | Original before | reuse | **missing** | both available | both available |
| 7CCO | Hybrid-La proposed after | reuse | reuse from13141 | both available | both available |
| 7CCN | Original before | reuse | **missing** | both available | both available |
| 7CCN | Whole-MACE-La step80 after | reuse | **missing** | both available | both available |

The three missing POINTs explicitly record `compute_force=false`; each linked raw MACE receipt has `forces_eV_A=null` and `force_convention=null`. This is an archived request choice, not unsupported Dy force functionality or a failed force calculation. The pinned `NativeMACE.evaluate` accepts physical Dy multiplicity6 and `compute_force=True`; the existing 7CCO after-Dy request already used that route successfully.

`INVENTORY.json` binds the exact two complete collectors:

- 7CCO job13424: RESULT `9229f7bebb533dd44a024bb80eca10eded2b3a64cc3972b245d65d82b5ecda97`, config `8f1446fd601f5dceb91a3f6bce290787a9fad10b89270124472bfaf02ddeeb81`, adapter `06e40d37a44c16446a062bbc479b6a49db70c72e6d16057a2fd7f7278916ee0f`.
- 7CCN job15173: RESULT `fff95c2f235947acb24e91991915b320731992f48c2b3d6c9ca3912a6e927f38`, config `a9ebe125029d581ecd18ed7f270ee5a625e73df547242770e90fd2c51427282a`, adapter `c0cca6dd4f44619fc2e3201f024b7913f6b7e9d8e8c4fca694fe17bd942d0ecf`.

Each of the eight inventory rows includes its admitted native row, full force report, physical IDs, geometry/task/XYZ pins, raw MACE receipt and available force artifacts, full reference state, native output/engrad/receipt, and saved solvent receipt/artifacts. The original 7CCO model completion13142 supplies both after MACE force sets and the supplemental after-Dy solvent; its original preflight solvent inventory is not rewritten.

The current and historical already-bound whole-pool menus exposed16 unique POINT pins. Exact source/metal/physical-hash matches for all three missing cells remain energy-only. This is a bounded cache check of those menus, not a claim to have searched every workspace or compared geometries numerically.

## Recipe and implementation constraints

Reuse original whole-force CONFIG `d3c28300c8c7ef88ef747b733603a344790b932db56aa53b0ccd4e999a7a4023` and native helper `0e48b8c6c2f271512f2c445f826036d8a68f7dca0aabc5923697ab3ecb7d9896`. They specify the existing checkpoint, float64, omol head, torch2.8.0/mace0.3.16/ase3.29.0 and **1e−7 eV** energy replay tolerance. That tolerance is reproducibility admission, not an accuracy criterion. Keep charge−2,256 atoms, target255, and physical multiplicities La1/Dy6.

The force-report labels differ between source adapters; explicit labels are recorded in each inventory row. Preserve original directions:7CCO full deformation chord and initial hybrid Phe1-CA force;7CCN full deformation chord only. No new direction is selected.

The reference solvent-gradient producer uses `(raw * HA) / BOHR`; whole-model force production uses `-raw * (HA / BOHR)`. Preserve each archived representation and operation order. Exact replay must not demand bitwise equality after reassociation. Native MACE forces are eV/Å; the frozen conversion is23.06054783061903 kcal/mol per eV.

## Scientific value and limits

The finite comparison can determine whether the existing energy agreement for Dy conceals force errors, and whether those errors cancel in the four matched Dy-minus-La force differences. This changes the evidence available for choosing a mechanics model; another La-only ranking would not answer it.

MACE's physical Dy6 conditioning and PBE0's effective-singlet frozen4f9 ECP55 representation remain distinct. La uses ECP46. Both reference models omit SOC; the vacuum dianion, conditional termini, finite-basis/state and printed SCF residual caveats persist. These conditional comparisons are meaningful under the already declared representations, but neither the three force calls nor the eight-cell coverage qualifies a physical ground state or validates the shared solvent. Native analytic-gradient admission remains separate from derivative qualification.

The original hybrid after-Dy result remains unavailable. The original7CCN step80 pair is retained; steps92/172/725 and live15545 are outside this inventory. No closed hybrid, electronic-base or Mex derivative branch is reopened. Natural-transfer mechanics/state gates remain unfulfilled.

## Review scope

Source and stdlib metadata/hash reads only:183 explicit non-checkpoint artifact pins verified, all eight reported reference force tables have256 ordered IDs, and raw request availability matches the original collectors. The large model checkpoint pin was copied from the qualified configuration and was not rehashed here; native runtime must verify it. No scientific imports, vector arithmetic, geometry recomputation, model loads, new calls, submissions or changes to original records occurred.

Inventory SHA256: `e2c196d1b5961cc74d50647e0565780928c3c1c103b282d563833f4b7e003957`.
